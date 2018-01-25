// Class for header information
class HeaderInfo{
public:
  //void clear(); // Reset header variables
  uint8_t id = 0; // FEM ID
  uint8_t slot = 0; // FEM slot in crate
  bool test = false; // Test mode flag
  bool overflow = false; // Overflow flag
  bool full = false; // Full flag
  uint32_t nwords = 0; // Number of words in frame
  uint32_t event = 0; // Event number
  uint32_t frame = 0; // Frame number
  uint32_t checksum = 0; // Checksum

  // Crosscheck header information
  uint32_t wordcount = 0; // Actual (counted) number of words in frame
  uint32_t mychecksum = 0; // Manual checksum

  // No members for the 6th pair of headers since it may be buggy/undefined yet
  void clear(){
    id = 0;
    slot = 0;
    test = false;
    overflow = false;
    full = false;
    nwords = 0;
    event = 0;
    frame = 0;
    checksum = 0;
    wordcount = 0;
    mychecksum = 0;
  };
};
