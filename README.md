# DecoderTools
Tools for decoding and analyzing the TPC data produced by the Nevis boards
The decoder will decode the binary file and split them into small root files, then the analyzer will analyze these root files, finally comparison.exe will compare the spectrum with reference run. (don't forget to hadd root files together before comparing)


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
To run analyzer on condor, we need to transfer files to /share/ disk.

where your_decoded_nevis_tpc_file.root corresponds to a run taken with BNL electronics in channel-map mode and your_bnl_pin_mapping.txt is a text file provided by BNL.
