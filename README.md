# DecoderTools
Tools for decoding and analyzing the TPC data produced by the Nevis boards

To compile the tools, run
```
./make.sh
```
To decode a binary file, run
```
./decoder.exe your_nevis_tpc_binary_file.dat
```
To plot a decoded file, run
```
./plotter.exe your_decoded_nevis_tpc_file.root
```
To create the channel map, run
```
./channel_mapper.exe your_decoded_nevis_tpc_file.root your_bnl_pin_mapping.txt
```
where your_decoded_nevis_tpc_file.root corresponds to a run taken with BNL electronics in channel-map mode and your_bnl_pin_mapping.txt is a text file provided by BNL.