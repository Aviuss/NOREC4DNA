import os
import json
import math
import bisect
import struct
import random
import argparse
# from pathos.multiprocessing import ProcessingPool as Pool
# from pathos import multiprocessing
# import dill
import progressbar
import multiprocessing
from functools import partial
from math import ceil, floor

from norec4dna.Encoder import Encoder
from norec4dna.LTEncoder import LTEncoder
from norec4dna.Packet import ParallelPacket
from norec4dna.RU10Encoder import RU10Encoder
from norec4dna.OnlineEncoder import OnlineEncoder
from norec4dna.rules.FastDNARules import FastDNARules
from norec4dna.helper.RepeatedTimer import RepeatedTimer
from helpful_scripts.automatedfindminimum import AutomatedFindMinimum
from norec4dna.distributions.OnlineDistribution import OnlineDistribution
from norec4dna.distributions.RaptorDistribution import RaptorDistribution
from norec4dna.ErrorCorrection import nocode, crc32, reed_solomon_encode, dna_reed_solomon_encode
from norec4dna.helper import should_drop_packet, split_file, number_to_base_str, find_ceil_power_of_four, \
    merge_folder_content
from norec4dna.distributions.ErlichZielinskiRobustSolitonDisribution import ErlichZielinskiRobustSolitonDistribution

DEFAULT_ID_LEN_FORMAT = "H"
DEFAULT_NUMBER_OF_CHUNKS_LEN_FORMAT = "I"
ONLINE_QUALITY = 7
ONLINE_EPS = 0.068

DEFAULT_CONFIG = {
    "filename": None,  # the file to encode (required)
    "repair_symbols": 2,  # number of repair symbols for ReedSolomon
    "crc_size": 1,  # size of the crc in bytes: 1=CRC-8, 2=CRC-16, 4=CRC-32, 8=CRC-64
    "list_size": 1000,  # total number of packets to generate, split across cores [ignored if sequential is set to true]
    "out_size": 1000,  # number of packets to save
    "chunk_size": 0,  # chunk size (0 = infer from number of chunks)
    "number_of_chunks": 300,  # number of chunks [ignored if chunk_size is set to value != 0]
    "sequential": False,
    "cores": None,  # number of worker processes (None = all cpu cores)
    "split_input": 1,  # number of sub-codes to split input file into
    "error_correction": "reedsolomon",  # nocode, crc or reedsolomon
    "insert_header": False,
    "method": "RU10",  # RU10, Online or LT
    "optimization": False,
    "overhead": 0.2,  # overhead to use for the optimization
    "overhead_factor": 0.0,  # factor to exceed the overhead with if necessary to optimize
    "error_prob_factor": 0.1,  # factor for the max error_prob for new packets based on the avg of the existing ones
    "plot": False,  # plot the results of the optimization
    "mode1bmp": False,  # input is an image, save as b/w bmp without header (applies to RU10 only)
    "drop_above": 1.0,  # instantly drop packet if error-prob above set threshold
    "seed_size_str": DEFAULT_ID_LEN_FORMAT,  # struct-string for seed: I,H,B (limits #packets in sequential mode)
    "xor_by_seed": False,
    "id_spacing": 0,
}

# DEFAULT_CHUNK_SIZE = 34
# NUMBER_OF_PACKETS_TO_CREATE = 655360

# global counter for progressbar
# counter = None
progress_bar = None


def run(seq_seed=None, file='logo.jpg', repair_symbols=2, insert_header=False,
        error_correction="reedsolomon", save_number_of_chunks_in_packet=False, l_size=1000, while_count=1000,
        chunk_size=0, number_of_chunks=300, prepend="", append="", seed_len_format=DEFAULT_ID_LEN_FORMAT,
        number_of_chunks_len_format=DEFAULT_NUMBER_OF_CHUNKS_LEN_FORMAT, method='RU10',
        mode1bmp=False, drop_above=0.4, packets_to_create=None, xor_by_seed=False, id_spacing=0, custom_dist=None,
        crc_size=1):
    # global counter
    if chunk_size != 0:
        number_of_chunks = Encoder.get_number_of_chunks_for_file_with_chunk_size(file, chunk_size)
    dna_rules = FastDNARules()
    if packets_to_create is None:
        packets_to_create = math.pow(2, 8 * struct.calcsize(seed_len_format))
    rules = dna_rules
    dist = get_dist(method, number_of_chunks)
    error_correction = instantiate_error_correction(error_correction, repair_symbols, crc_size)
    
    if custom_dist is not None:
        dist.f = custom_dist
        dist.d = [x for x in range(0, 41)]

    if method == 'RU10':
        x = RU10Encoder(
            file, number_of_chunks, dist,
            chunk_size=chunk_size,
            insert_header=insert_header,
            rules=rules,
            error_correction=error_correction,
            id_len_format=seed_len_format,
            number_of_chunks_len_format=number_of_chunks_len_format,
            save_number_of_chunks_in_packet=save_number_of_chunks_in_packet,
            mode_1_bmp=mode1bmp,
            prepend=prepend,
            append=append,
            xor_by_seed=xor_by_seed,
            id_spacing=id_spacing
        )
        x.prepare()
    elif method == 'LT':
        x = LTEncoder(file, number_of_chunks, dist,
            chunk_size=chunk_size,
            insert_header=insert_header,
            rules=rules,
            error_correction=error_correction, number_of_chunks_len_format=number_of_chunks_len_format,
            id_len_format=seed_len_format, save_number_of_chunks_in_packet=save_number_of_chunks_in_packet
        )
        x.prepareEncoder()
    elif method == 'Online':
        number_of_chunks = dist.get_size()
        x = OnlineEncoder(file, number_of_chunks, dist, ONLINE_EPS, ONLINE_QUALITY,
            error_correction=error_correction,
            quality_len_format="B",
            insert_header=False,
            check_block_number_len_format=seed_len_format,
            number_of_chunks_len_format=number_of_chunks_len_format,
            rules=rules,
            save_number_of_chunks_in_packet=False
        )
        x.prepare()
    else:
        raise NotImplementedError("Choose: RU10, LT or Online")

    i = 0
    tmp_list = []
    while i < while_count:
        if seq_seed is not None:
            if seq_seed + i >= packets_to_create:
                break
            packet = x.create_new_packet(seed=seq_seed + i)
        else:
            packet = x.create_new_packet()

        _ = should_drop_packet(rules, packet)
        if packet.error_prob <= drop_above and (len(tmp_list) < l_size or packet.error_prob < tmp_list[-1].error_prob):
            if packet not in tmp_list:
                bisect.insort_left(tmp_list, packet)
            else:
                elem = next((x for x in tmp_list if x == packet), None)
                if packet < elem:
                    tmp_list.remove(elem)
                    del elem
                    bisect.insort_left(tmp_list, packet)
            if len(tmp_list) > l_size:
                for ele1m in tmp_list[l_size + 1:]:
                    del ele1m
                tmp_list = tmp_list[:l_size]

        else:
            del packet
        i += 1
        # += operation is not atomic, so we need to get a lock:
        # with counter.get_lock():
        #    counter.value += 1
    # save_packets_fasta(tmp_list, out_file=method + "_out_partial", file_ending="." + method + "_DNA",
    #                   clear_output=False)
    
    conf = {'error_correction': error_correction, 'repair_symbols': repair_symbols,
            # 'number_of_splits': _number_of_splits,
            'find_minimum_mode': True, 'seq_seed': seq_seed}
    
    # x.save_config_file(conf, section_name=method + "_" + file)
    if x.progress_bar is not None:
        x.progress_bar.finish()
    return [ParallelPacket.from_packet(p) for p in tmp_list]

def run_worker(seq_seed, while_count, **kwargs):
    return run(seq_seed=seq_seed, while_count=while_count, **kwargs)

def save_packets_fasta(packets, out_file, file_ending, seed_is_filename=True):
    i = 0
    if not os.path.exists(out_file):
        os.makedirs(out_file)
    out_str = out_file + "/" + str(random.randint(1, 80000000000)) + "_" + str(
        random.randint(1, 800000000000)) + ".fasta"
    with open(out_str, "w") as f:
        for packet in packets:  # sorted(packets, key=lambda elem: (elem.error_prob, elem.__hash__())):
            if seed_is_filename:
                i = packet.id
            e_prob = (str(ceil(packet.error_prob * 100)) + "_") if packet.error_prob is not None else ""
            f.write(">" + e_prob + str(i) + file_ending + "\n" + packet.get_dna_struct(True) + "\n")
            i += 1
    print(f"Saved result at: %s" % out_str)


def reduce_lists(base, input_list=None, l_size=100):
    if input_list is None and len(base) == 2:
        input_list = base[1]
        base = base[0]
    base = base[:l_size]
    for packet in input_list:
        if packet.error_prob < base[-1].error_prob:
            if packet not in base:
                bisect.insort_left(base, packet)
            else:
                elem = next((x for x in base if x == packet), None)
                if packet < elem:
                    base.remove(elem)
                    bisect.insort_left(base, packet)
            if len(base) > l_size:
                base = base[:l_size]
        else:
            # we can skip the rest of the current result-row since all remaining packets have a higher error_prob
            # than the worst packet of the current answer
            return base
    return base


# def init_mp(args):
#    global counter
#    counter = args


def create_progress_bar(max_value):
    widgets = [
        progressbar.Percentage(),
        progressbar.Bar(),
        ' Correct: ',
        progressbar.Counter(),
        ', ',
        progressbar.AdaptiveETA(), ' ',
        progressbar.Timer()
    ]
    return progressbar.ProgressBar(max_value=max_value, widgets=widgets, max_error=False,
                                   redirect_stdout=True).start()


def update_progressbar():
    pass
    # progress_bar.update(counter.value)


# multiprocess and merge all created lists at the end
def main(filename="logo.jpg", repair_symbols=2, while_count=1000, out_size=1000, chunk_size=0, number_of_chunks=300,
         sequential=False, prepend="", append="", insert_header=False,
         seed_len_format=DEFAULT_ID_LEN_FORMAT,
         method='RU10', mode1bmp=False, drop_above=0.4,
         packets_to_create=None, xor_by_seed=False, id_spacing=0, custom_dist=None, error_correction="reedsolomon",
         crc_size=1, cores=None):
    # global progress_bar, counter
    if packets_to_create is None:
        packets_to_create = math.pow(2, 8 * struct.calcsize(seed_len_format))
    if cores is None:
        cores = multiprocessing.cpu_count()
    # counter = multiprocessing.Value('i', 0)
    p = multiprocessing.Pool(cores)  # , initializer=init_mp, initargs=(counter,))
    if sequential:
        stepsize = packets_to_create / cores
        param = [(int(floor(i * stepsize)), int(ceil(stepsize)) + 1) for i in range(cores)]
    else:
        param = [(None, while_count // cores + (1 if i < while_count % cores else 0)) for i in range(cores)]

    progress_bar = create_progress_bar(packets_to_create)
    progress_bar.update(0)
    rt = RepeatedTimer(2, update_progressbar, )
    a = p.starmap(
        partial(run_worker, file=filename, repair_symbols=repair_symbols, l_size=out_size,
                chunk_size=chunk_size, number_of_chunks=number_of_chunks, prepend=prepend, append=append,
                insert_header=insert_header, seed_len_format=_seed_size_str, method=method, mode1bmp=mode1bmp,
                drop_above=drop_above, xor_by_seed=xor_by_seed, id_spacing=id_spacing, custom_dist=custom_dist,
                error_correction=error_correction, crc_size=crc_size), param)
    rt.stop()
    progress_bar.finish()
    print("Merging results...")
    progress_bar.start(len(a))
    progress_bar.update(0)
    while len(a) > 1:
        if len(a) % 2 != 0:
            a.append([])
        tmp = []
        for i in range(0, len(a) - 1, 2):
            tmp.append([a[i], a[i + 1]])
        a = p.map(partial(reduce_lists, input_list=None, l_size=out_size), tmp)
        progress_bar.update(progress_bar.max_value - len(a))
    a = a[0]
    rt.stop()

    out_file = "output"
    file_ending = "_joined" + ("_" + append if append else "") + "." + method + "_DNA"
    save_packets_fasta(a, out_file=out_file, file_ending=file_ending)
    progress_bar.finish()
    return a, out_file, file_ending


def get_dist(_method, _number_of_chunks):
    if _method == 'RU10':
        dist = RaptorDistribution(_number_of_chunks)
    elif _method == 'LT':
        dist = ErlichZielinskiRobustSolitonDistribution(_number_of_chunks, seed=2)
    elif _method == 'Online':
        dist = OnlineDistribution(ONLINE_EPS)
    else:
        raise NotImplementedError("Choose: RU10, LT or Online")
    return dist


def init_optimization(packets, error_correction, dist):
    return AutomatedFindMinimum(packets, error_correction=error_correction, dist=dist)


def instantiate_error_correction(err_corr_name, no_repair_symbols, crc_size):
    if err_corr_name == "nocode":
        return nocode
    elif err_corr_name == "crc":
        crc_len_format = {1: "B", 2: "H", 4: "L", 8: "Q"}[crc_size]
        return lambda x: crc32(x, crc_len_format)
    elif err_corr_name == "reedsolomon":
        return lambda x: reed_solomon_encode(x, no_repair_symbols)
    else:
        raise RuntimeError(
            "Selected Error Correction not supported, choose: 'nocode', 'crc', 'reedsolomon'")


def load_config(config_path):
    with open(config_path) as f:
        config = json.load(f)
    if not isinstance(config, dict):
        raise ValueError("Config file must contain a json object")
    unknown_keys = set(config) - set(DEFAULT_CONFIG)
    if unknown_keys:
        raise ValueError("Unknown config parameter(s): " + ", ".join(sorted(unknown_keys)))
    config = {**DEFAULT_CONFIG, **config}
    if config["filename"] is None:
        raise ValueError("Config parameter 'filename' (the file to encode) is required")
    if config["error_correction"] not in ("nocode", "crc", "reedsolomon"):
        raise ValueError("Config parameter 'error_correction' must be one of: nocode, crc, reedsolomon")
    if config["crc_size"] not in (1, 2, 4, 8):
        raise ValueError("Config parameter 'crc_size' must be one of: 1, 2, 4, 8")
    if config["cores"] is not None and (not isinstance(config["cores"], int) or isinstance(config["cores"], bool)
                                        or config["cores"] < 1):
        raise ValueError("Config parameter 'cores' must be a positive integer or null")
    return config


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("config", metavar="config", type=str,
                        help="path to a json file with the parameters (see DEFAULT_CONFIG for keys and defaults)")
                        
    config = load_config(parser.parse_args().config)
    _input_file = config["filename"]
    _repair_symbols = config["repair_symbols"]
    _list_size = config["list_size"]
    _out_size = config["out_size"]
    _chunk_size = config["chunk_size"]
    _number_of_hunks = config["number_of_chunks"]
    _sequential = config["sequential"]
    _cores = config["cores"]
    _number_of_splits = config["split_input"]
    _insert_header = config["insert_header"]
    _error_correction = config["error_correction"]
    _drop_above = config["drop_above"]
    _seed_size_str = config["seed_size_str"]
    _method = config["method"]
    _mode1bmp = config["mode1bmp"]

    _optimization = config["optimization"]
    _overhead = config["overhead"]
    _overhead_factor = config["overhead_factor"]
    _error_prob_factor = config["error_prob_factor"]
    _plot = config["plot"]
    _xor_by_seed = config["xor_by_seed"]
    _id_spacing = config["id_spacing"]
    _crc_size = config["crc_size"]

    if _chunk_size != 0:
        _number_of_hunks = Encoder.get_number_of_chunks_for_file_with_chunk_size(_input_file, _chunk_size)
    if _list_size < _out_size:
        print("Warning: A list size smaller than the out_size may result in non-optimal results.")
    if _optimization:
        print("Warning: Please make sure to generate enough packets for the automated optimization.")

    if not _optimization:
        if _number_of_splits > 1:
            input_files = split_file(_input_file, _number_of_splits)
            power_of_four = find_ceil_power_of_four(len(input_files))
            print("Spltting input into {} sub files. We need to prepend/append {} base(s)".format(len(input_files),
                                                                                                  power_of_four))
            prepend_matching = {input_files[i]: number_to_base_str(i, power_of_four) for i in range(len(input_files))}
        else:
            input_files = [_input_file]
            prepend_matching = {_input_file: ""}
        for _input_file in input_files:
            try:
                print("File to encode: " + str(_input_file))
                # def callab():
                main(_input_file, _repair_symbols, _list_size, _out_size, _chunk_size, _number_of_hunks,
                     _sequential, method=_method, append=prepend_matching[_input_file],
                     insert_header=_insert_header, seed_len_format=_seed_size_str, drop_above=_drop_above,
                     xor_by_seed=_xor_by_seed, id_spacing=_id_spacing, error_correction=_error_correction,
                     crc_size=_crc_size, cores=_cores)

                # print(timeit.repeat(callab, number=1111, repeat=5))
            except Exception as ex:
                print(ex)
                raise ex
                # list_fds()
        # if len(input_files) > 1:
        #    merge_folder_content("split_" + os.path.basename(_input_file), _input_file + "_combined_split_output",
        #                         append_folder_name=True,
        #                         clear_dest_folder=True)
    else:
        par_packets, _out_file, _file_ending = main(_input_file, _repair_symbols, _list_size, _out_size, _chunk_size,
                                                    _number_of_hunks, _sequential, method=_method,
                                                    seed_len_format=_seed_size_str, drop_above=_drop_above,
                                                    error_correction=_error_correction, crc_size=_crc_size,
                                                    cores=_cores)
        distribution = get_dist(_method, _number_of_hunks)
        optimizer = init_optimization(par_packets,
                                      instantiate_error_correction(_error_correction, _repair_symbols, _crc_size),
                                      distribution)
        optimized_packets = optimizer.automated_optimization(_overhead, _overhead_factor, _error_prob_factor,
                                                             plot=_plot)
        save_packets_fasta(optimized_packets, out_file=_out_file + "optimized/", file_ending=_file_ending)
