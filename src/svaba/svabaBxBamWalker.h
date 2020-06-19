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

class svabaBxBamWalker : private SeqLib::BamReader {
  /* Reads a BAM file that was produced by the lariat aligner. This file must be
   prepared by flipping the chromosome and BX tag with bxtools convert. Then
   this file must be sorted and indexed by the BX tag using samtools. This
   allows fast retrieval of reads having a specific bx tag.
  */

public:
  /* bx_bam_path: BAM file indexed by bx tag. */
  svabaBxBamWalker(const std::string &bx_bam_path);

  /*  */
  svabaReadVector fetchReadsByBxBarcode(const BxBarcode &bx_barcode);
  svabaReadVector fetchReadsByBxBarcode(const std::set<BxBarcode> &bx_barcodes);
  SeqLib::BamHeader Header() const;

  static std::set<BxBarcode> collectBxBarcodes(const svabaReadVector &reads);
};

#endif
