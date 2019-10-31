#! /bin/bash

cd share/LHAPDF/

LHAPDF_version=6.1.6
PDFs="NNPDF30_lo_as_0118 NNPDF30_nlo_as_0118 NNPDF30_nnlo_as_0118"

for PDF in $PDFs; do
  echo "Installing $PDF PDF set"
  curl -O -L https://www.hepforge.org/archive/lhapdf/pdfsets/$LHAPDF_version/$PDF.tar.gz
  tar xf $PDF.tar.gz
  rm $PDF.tar.gz
done;


