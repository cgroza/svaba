#include "svabaBxBamWalker.h"
#include <iostream>
#include <stdexcept>
#include <system_error>

svabaBxBamWalker::svabaBxBamWalker(const std::string &bx_bam_path,
                                   const std::string _prefix,
                                   bool _weird_reads_only)
    : prefix(_prefix), weird_reads_only(_weird_reads_only) {

  BamReader();
  Open(bx_bam_path);
}

svabaBxBamWalker::svabaBxBamWalker() : weird_reads_only(true) { BamReader(); }

svabaReadVector
svabaBxBamWalker::fetchReadsByBxBarcode(const BxBarcode &bx_barcode) {
  svabaReadVector read_vector;
  // We must convert the stirng barcode into an index ID. We retrieve this ID
  // from the BAM header.
  SeqLib::BamHeader header = Header();
  SeqLib::GenomicRegion bx_region(bx_barcode, "1", "2", header);

  bool success_seek = SetRegion(bx_region);
  if (!success_seek) {
    auto ex =
        new std::invalid_argument("Could not set region to " + bx_barcode);
    std::cerr << ex->what() << std::endl;
    throw ex;
  }

  // BX tags are sorted and indexed in large contiguous blocks within the BAM.
  // Each BX tag is a region. Therefore, we stop when we exhaust the region.
  while (true) {
    // Careful. We cannot reuse BamRecords since we must avoid pushing shallow
    // copies into the svabaReadVector.
    SeqLib::BamRecord bx_record;

    // Are we still within the same BX block?
    if (GetNextRecord(bx_record)) {
      if (isBxReadWeird(bx_record))
        read_vector.push_back(svabaRead(bx_record, prefix));
    } else
      break;
  }
  return read_vector;
}

svabaReadVector svabaBxBamWalker::fetchReadsByBxBarcode(
    const std::set<BxBarcode> &bx_barcodes) {
  svabaReadVector all_reads;
  for (const BxBarcode &barcode : bx_barcodes) {
    std::string barcode_copy = barcode;
    // "-" gets translated to "_" during bxtools convert. We must correct this.
    std::replace(barcode_copy.begin(), barcode_copy.end(), '-', '_');
    svabaReadVector reads = fetchReadsByBxBarcode(barcode_copy);
    // move the reads into all reads
    std::move(reads.begin(), reads.end(), std::back_inserter(all_reads));
  }
  return all_reads;
}

std::set<BxBarcode>
svabaBxBamWalker::collectBxBarcodes(const svabaReadVector &reads) {
  std::set<BxBarcode> barcodes;
  for (const svabaRead &r : reads) {
    std::string bx_tag;
    // tag may not always be present
    if (r.GetZTag("BX", bx_tag)) {
      barcodes.insert(bx_tag);
    }
  }
  return barcodes;
}

bool svabaBxBamWalker::isBxReadWeird(SeqLib::BamRecord &r) {
  // look for unmapped reads, unpaired reads and poor alignments
  // NOTE: the ProperPair flag is part of the SAM specification and is set by
  // the aligner Includes information on insert size, mate read orientation
  // etc...
    bool isWeird = !r.PairMappedFlag() || !r.MappedFlag() ||
        (r.MapQuality() <= POOR_ALIGNMENT_MAX_MAPQ); // ||
        // !r.ProperPair();

    // report which flag was set
#ifdef DEBUG_SVABA_BX_BAM_WALKER
  if (isWeird)
    std::cerr << !r.PairMappedFlag() << !r.MappedFlag() << !r.ProperPair()
              << (r.MapQuality() < POOR_ALIGNMENT_MAX_MAPQ) << std::endl;
#endif
  return isWeird;
}
