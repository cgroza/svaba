#ifndef SVABA_BX_BAM_WALKER_H
#define SVABA_BX_BAM_WALKER_H

#include "SeqLib/BFC.h"
#include "SeqLib/BamReader.h"
#include "SeqLib/BamRecord.h"
#include "SeqLib/GenomicRegion.h"
#include "svabaRead.h"
#include <algorithm>
#include <set>
#include <string>
#include <vector>
#include <iostream>
#include <exception>

/* Let's distinguish regular strings from BxBarcodes in the source code */
typedef std::string BxBarcode;

class svabaBxBamWalker : public SeqLib::BamReader {
  /* Reads a BAM file that was produced by the lariat aligner. This file must be
   prepared by flipping the chromosome and BX tag with bxtools convert. Then
   this file must be sorted and indexed by the BX tag using samtools. This
   allows fast retrieval of reads having a specific bx tag.
  */

public:
  /* bx_bam_path: BAM file indexed by bx tag. */
    svabaBxBamWalker();
    svabaBxBamWalker(const std::string &bx_bam_path,
                     const std::string _prefix = "0000",
                     bool _weird_reads_only = true);

  /*  */
  svabaReadVector fetchReadsByBxBarcode(const BxBarcode &bx_barcode);
  svabaReadVector fetchReadsByBxBarcode(const std::set<BxBarcode> &bx_barcodes);

  static std::set<BxBarcode> collectBxBarcodes(const svabaReadVector &reads);
  std::string prefix;

  static bool isBxReadWeird(SeqLib::BamRecord &r);
  static const int POOR_ALIGNMENT_MAX_MAPQ = 30;

private:
    bool weird_reads_only;

};

#endif
