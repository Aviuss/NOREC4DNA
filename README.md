# Fork

this is a fork

# Encoding in LT

```
sudo docker run --rm -v "$PWD":/data:Z -v "$PWD/output":/norec4dna/output:Z norec4dna /data/find_minimum_config.json
```

```
docker run --rm \
-v ./[input_folder]:/data \
norec4dna [file to encode] \
--method LT \
```



# Other

`docker build . --tag norec4dna`

packet bytes = seed bytes + chunk_size + repair_symbols
sequence (nt) = 4 × packet bytes (1 byte = 4 bases)



```
docker run --rm \
-v "$PWD":/data:Z \
norec4dna README.md \
--method RU10 --chunk_size 40 --insert_header   --repair_symbols 2 --seed_size_str I --list_size 2000 --out_size 100
```



# fragments of old NOREC4DNA docs

## Find minimum

### Building

Build the docker container

### Parameters

First enter the filename of the file to generate the packets for.

`FILE (--parameters)`

The following parameters can be set:

`--repair_symbols=[no_symbols]`

The number of repair_symbols for ReedSolomon (default=2). This does only apply if --error_correction is set to reedsolomon

`--list_size=[size]`

Size of operational list per thread, inferred by the number cores if sequential is set to true. The list size should
always be greater than the out_size to ensure optimal results (default=1000).

`--out_size=[size]`

Number of packets to save after combining the lists and sorting them by the packets error_prob (default=1000).

`--chunk_size=[size]`

Size of chunks to split the file into, inferred from number of chunks and the filesize if not set (default=0).

`--sequential`

If set, all seed will be generated in a sequential matter. (Recommended!)

`--spare1core`

If activated, one core is not used for the calculation of the lists.

`--method=[RU10/Online/LT]`

Sets the method to generate the packets with. Available are RU10, Online and LT.

`--seed_size_str=[I,H,...]`

Set the struct-string for the seed field. See [https://docs.python.org/3/library/struct.html#format-characters](https://docs.python.org/3/library/struct.html#format-characters) for more information

`--drop_above`

Sets an upper-limit for the error probability. WARNING: This might reduce the total number of sequences returned!

## With optimization

`--optimization`

Activates the automated optimization of the chunk distribution in the packets with different options.

`--overhead=[overhead (0.1=10%)]`

Overhead to use for the optimization, where 0.1 means 10% additional packets based on the number of packets needed to
decode the file (default=0.1).

`--overhead_factor=[factor (0.1=10%)]`

If the overhead is not enough to optimize the packets, the overhead factor is a factor that allows exceeding the given
overhead to try to optimize the chunk distribution (default=0.0).

`--errorprob_factor=[factor (0.1=10%)]`

A factor for the maximum allowed error_prob of the additional packets based on the average error_prob of the packets
needed to decode (default=0.1). If set to 0.0 no more packets may be added since the packets with the lowest error_probs
were already used to decode the file.

`--plot`
Generates and saves different plots to show the results.

---

## Tools

### demo_*.py

Demo applications for fast en- and decoding of sequences.

### ru10_find_minimum_packets.py (Deprecated)
`--error_correction [nocode, crc, reedsolomon]`

Defines the error detection / correction algorithm to use per packet. (Default: nocode = no error-detection/correction)

`--split_input`

Sets the number of pre-splits to perform Default: 1 (= do not split the input file into multiple NOREC rounds)
WARNING: If set, this value _should_ be known during decoding (thus using a bruteforce approach this value might be reconstructed) 

`--store_as_fasta`

If set, stores the result in a .fasta file instead of one file per sequence

`--insert_header`

If set, besides the created chunks an additional header chunk will be added. This chunk stores the filename and the correct padding for the last chunk.
(Recommended!) WARNING: If not set, the reconstructed file will most likely be longer due to the \00-padding at the end.
 

### ConfigWorker.py
Allows easy en- and decoding used .ini files.
Since the supplied encoder can create such .ini files, this is especially useful for easy decoding.


### helpful_scripts

there are various more or less useful scripts inside `helpful_scripts/`


---


## Example
#### To try out NOREC4DNA you can use the demo_\*\_encode.py python scripts:

`python demo_raptor_encode.py Dorn --error_correction=reedsolomon --repair_symbols=3 --as_dna --insert_header`

###### this should create a new folder "RU10_Dorn" as well as an Dorn\_\*.ini file.

#### To decode the file from DNA one could either use demo\_\*\_decode.py:

`python demo_raptor_decode.py RU10_Dorn --use_header_chunk --error_correction=reedsolomon --repair_symbols=3 --number_of_chunks=145 (number as seen in the ini, unless --save_number_of_chunks was defined during encoding)`

#### or use the ConfigWorker.py:

`python ConfigWorker.py <name of the .ini-file>`

The decoded file will be saved as DEC_RU10_Dorn if no header-chunk was added during encoding, otherwise the file will be saved under the correct filename.
###### if the header-chunk was NOT used, the created file will have padding \00-bytes at the end. 
