#include "artdaq-core-mu2e/Data/CalorimeterTestDataDecoder.hh"

#include "TRACE/tracemf.h"

#include <algorithm>
#include <cmath>


namespace mu2e {

  CalorimeterTestDataDecoder::CalorimeterTestDataDecoder(DTCLib::DTC_SubEvent const& evt)
    : DTCDataDecoder(evt)
  {
    if (block_count() > 0)
      {
        auto dataPtr = dataAtBlockIndex(0);// Return pointer to beginning of DataBlock at given DataBlock index( returns type DTCLib::DTC_DataBlock)
        auto hdr = dataPtr->GetHeader(); // get the header
        // check that the subsystem is the calo and the version is correct:
        if (hdr->GetSubsystem() != DTCLib::DTC_Subsystem_Calorimeter || hdr->GetVersion() > 1)
          {
            TLOG(TLVL_ERROR) << "CalorimeterTestDataDecoder CONSTRUCTOR: First block has unexpected type/version " << hdr->GetSubsystem() << "/" << static_cast<int>(hdr->GetVersion()) << " (expected " << static_cast<int>(DTCLib::DTC_Subsystem_Calorimeter) << "/[0,1])";
          }
      }
  }

  CalorimeterTestDataDecoder::CalorimeterTestDataDecoder(std::vector<uint8_t> data)
    : DTCDataDecoder(data)
  {
    //event_.GetDataBlockCount() > 0
    if (block_count() > 0){
      auto dataPtr = dataAtBlockIndex(0);
      auto hdr = dataPtr->GetHeader();
      if (hdr->GetSubsystem() != DTCLib::DTC_Subsystem_Calorimeter || hdr->GetVersion() > 1){
        TLOG(TLVL_ERROR) << "CalorimeterTestDataDecoder CONSTRUCTOR: First block has unexpected type/version " << hdr->GetSubsystem() << "/" << static_cast<int>(hdr->GetVersion()) << " (expected " << static_cast<int>(DTCLib::DTC_Subsystem_Calorimeter) << "/[0,1])";
      }
    }
  }

  // Get Calo Hit Data Packet
  std::vector<std::pair<mu2e::CalorimeterTestDataDecoder::CalorimeterHitTestDataPacket, std::vector<uint16_t>>>* mu2e::CalorimeterTestDataDecoder::GetCalorimeterHitData(size_t blockIndex) const
  {
    std::vector<std::pair<mu2e::CalorimeterTestDataDecoder::CalorimeterHitTestDataPacket, std::vector<uint16_t>>> *output = new std::vector<std::pair<mu2e::CalorimeterTestDataDecoder::CalorimeterHitTestDataPacket, std::vector<uint16_t>>>();
  
    // get data block at given index
    DTCLib::DTC_DataBlock const * dataBlock = dataAtBlockIndex(blockIndex);
    if (dataBlock == nullptr) return output;

    DTCLib::DTC_DataHeaderPacket* blockHeader = dataBlock->GetHeader().get();
    size_t blockSize = dataBlock->byteSize;
    size_t nPackets = blockHeader->GetPacketCount();
    size_t dataSize = blockSize - 16;

    auto blockDataPtr = dataBlock->GetData();

    if (nPackets == 0){ //Empty packet
      TLOG(TLVL_WARNING) << "CalorimeterTestDataDecoder::GetCalorimeterHitData : no packets -- disabled ROC?\n";
      return output;
    }
    
    auto blockPos = reinterpret_cast<const uint8_t*>(blockDataPtr); //byte position in block (multiple of 16)
	  while(blockPos < reinterpret_cast<const uint8_t*>(blockDataPtr) + dataSize){ //until the end of this block      

      mu2e::CalorimeterTestDataDecoder::Data12bitReader reader(reinterpret_cast<const uint16_t*>(blockPos));
      
      //Make sure first word is 0xAAA
      if (reader[0] != 0xAAA){
        TLOG(TLVL_ERROR) << "CalorimeterTestDataDecoder::GetCalorimeterHitData : BeginMarker is " << std::hex << reader[0] << std::dec << " instead of 0xAAA\n";
        return output;
      }
      
      //std::cout<<"Searching for 0xFFF...\n";
      uint nWordsMax = ((dataSize * 8) / 12);
      int lastSampleMarkerIndex = -1;
      for (uint i=4; i<nWordsMax; i++){ //waveform starts from 5th word
        if (reader[i] == 0xFFF){
          //std::cout<<"FOUND 0xFFF at index "<<i<<"\n";
          lastSampleMarkerIndex = i;
          break;
        }
      }

      //Check if last marker was found
      if (lastSampleMarkerIndex == -1){
        TLOG(TLVL_ERROR) << "CalorimeterTestDataDecoder::GetCalorimeterHitData : LastSampleMarker 0xFFF not found in the payload!" << std::endl;
        return output;
      }

      //std::cout<<"Dumping data in 12-bit words:\n";
      //size_t nWordsTot = lastSampleMarkerIndex+6;
      //for (uint i=0; i<nWordsTot; i++){
      //  std::cout<< std::hex << std::setw(3) << std::setfill('0') << reader[i] << " " << std::dec;
      //  if (i==3 || i==nWordsTot-6) std::cout<<"\n";
      //}
      //std::cout<<"\n";

      //Create output
	  	output->emplace_back(mu2e::CalorimeterTestDataDecoder::CalorimeterHitTestDataPacket(), std::vector<uint16_t>());

      //Before waveform
      output->back().first.BeginMarker = reader[0];
      output->back().first.BoardID = reader[1];
      output->back().first.ChannelID = reader[2];
      output->back().first.InPayloadEventWindowTag = reader[3];

      //waveform
      size_t nSamples = lastSampleMarkerIndex-4;
      //std::cout<<"Saving waveform with "<<nSamples<<" samples\n";
		  output->back().second.resize(nSamples);
      for (uint i=0; i<nSamples; i++){
        output->back().second[i] = reader[4+i];
      }

      //After waveform
      output->back().first.LastSampleMarker = reader[lastSampleMarkerIndex];
      output->back().first.ErrorFlags = reader[lastSampleMarkerIndex+1];
      output->back().first.Time = (reader[lastSampleMarkerIndex+2] << 12) | reader[lastSampleMarkerIndex+3] ;
      output->back().first.IndexOfMaxDigitizerSample = reader[lastSampleMarkerIndex+4];
      output->back().first.NumberOfSamples = reader[lastSampleMarkerIndex+5];    

      //Advance to the next 16-byte packet
      float hitByteSize = nSamples*1.5 + sizeof(output->back().first);
      uint8_t hitPackets = uint8_t(std::ceil(hitByteSize/16)); //number of 16-byte packets this hit occupied
      blockPos += hitPackets*16; //advance by 16 bytes per packet
    }

    return output;
  }


  // Get Calo Counters Data Packet
  /*
  std::vector<std::pair<mu2e::CalorimeterTestDataDecoder::CalorimeterCountersDataPacket, std::vector<uint32_t>>>* mu2e::CalorimeterTestDataDecoder::GetCalorimeterCountersData(size_t blockIndex) const
  {
    // a pair is created mapping the data packet to set of hits (?)
    std::vector<std::pair<mu2e::CalorimeterTestDataDecoder::CalorimeterCountersDataPacket, std::vector<uint32_t>>> *output = new std::vector<std::pair<mu2e::CalorimeterTestDataDecoder::CalorimeterCountersDataPacket, std::vector<uint32_t>>>();
  
    // get data block at given index
    auto dataPtr = dataAtBlockIndex(blockIndex);
    if (dataPtr == nullptr) return output;

    // check size of hit data packet
    static_assert(sizeof(mu2e::CalorimeterTestDataDecoder::CalorimeterCountersDataPacket) % 2 == 0);
    
    //Get data in 32-bit words
    size_t headerSize = sizeof(*event_.GetHeader());
    size_t payloadSize = data_.size() - headerSize;
    auto payload32bitPtr = reinterpret_cast<uint32_t const*>(dataPtr + headerSize);
    std::vector<uint32_t> payload32bitData;
    for (uint word=0; word<payloadSize/4; word++){
      payload32bitData.push_back(*(payload32bitPtr + word));
    }

    //uint nSamples = payload32bitData.size();

    CalorimeterCountersDataPacket countersPacket(payload32bitData);

    output->emplace_back(countersPacket, payload32bitData);

    return output;
  }
  */

} // namespace mu2e
