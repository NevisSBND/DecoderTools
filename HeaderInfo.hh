// Class for header information
class HeaderInfo{
public:
  //void clear(); // Reset header variables
  size_t id = 0; // FEM ID
  size_t slot = 0; // FEM slot in crate
  uint32_t nwords = 0; // Number of words in frame
  uint32_t event = 0; // Event number
  uint32_t frame = 0; // Frame number
  uint32_t checksum = 0; // Checksum
  // No members for the 6th pair of headers since it may be buggy/undefined yet
  void clear(){
    id = 0;
    slot = 0;
    nwords = 0;
    event = 0;
    frame = 0;
    checksum = 0;
  };
};
