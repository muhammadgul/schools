# clone the lyon setup
git clone https://github.com/IPNL-CMS/HTTMadgraphDocumentation.git madgraph
cd madgraph

# install the latest PDF
curl -O -L https://www.hepforge.org/archive/lhapdf/LHAPDF-6.1.6.tar.gz
tar xf LHAPDF-6.1.6.tar.gz
rm LHAPDF-6.1.6.tar.gz
cd LHAPDF-6.1.6
./configure --prefix=$PWD/..
make -j8
make install
cd ..
#in case of lhapdf path problem use this export PATH=$PWD/bin:$PATH
# or this export LD_LIBRARY_PATH=$PWD/lib:$LD_LIBRARY_PATH
#download the proper pdf sets
./install_PDF_set.sh
#copy madgraph recent version
wget http://madgraph.physics.illinois.edu/Downloads/MG5_aMC_v2.6.0.tar.gz
tar xf MG5_aMC_v2.6.0.tar.gz
rm MG5_aMC_v2.6.0.tar.gz
#set the pdf path
path=$PWD/bin/lhapdf-config
sed "s@# lhapdf = lhapdf-config@$path@" MG5_aMC_v2_6_0/input/mg5_configuration.txt > MG5_aMC_v2_6_0/input/newmg5_configuration.txt
mv MG5_aMC_v2_6_0/input/newmg5_configuration.txt MG5_aMC_v2_6_0/input/mg5_configuration.txt
#download the SysCalc
wget http://madgraph.physics.illinois.edu/Downloads/SysCalc_V1.1.7.tar.gz
tar xf SysCalc_V1.1.7.tar.gz
rm SysCalc_V1.1.7.tar.gz
#cd SysCalc/src
#sed "s@lhapdf-config@$path@" Makefile > newMakefile && mv newMakefile Makefile
#make

# copy model to madgraph dir
cp -r model/Massive_Higgs_UFO MG5_aMC_v2_6_0/models/

#Generate template
#You can now generate your template using Madgraph. Let try to generate a scalar Higgs with a mass of 500 GeV, a width of 50 GeV and a coupling to ttbar of 1. Here are the different syntax to generate signal and interference separately.

cd MG5_aMC_v2_6_0/
./bin/mg5_aMC

#Generate signal only

MG5_aMC> import model Massive_Higgs_UFO
MG5_aMC> generate p p > h0 > t t~
MG5_aMC > output template_pp_h0_tt_S

#Generate interference only

MG5_aMC> import model Massive_Higgs_UFO
MG5_aMC> generate p p > t t~ / a0 HIG=2 HIG^2==2
MG5_aMC > output template_pp_h0_tt_i

# Generate Signal+Interferene

MG5_aMC> import model Massive_Higgs_UFO
MG5_aMC> generate g g > t t~ / a0 HIG=2 
MG5_aMC > output template_gg_h0_tt_S_i

#WARNING! Using this last squared order syntax, it is not possible to then decay the particles, neither in Madgraph nor in MadSpin, as reported in Madgraph launchpad. As you can see in the bug report, there is no "elegant" solution to this problem, and we have found dirty hacks to get around this problem which is the following: Generate Signal+Interferene template using generate g g > t t~ / a0 HIG=2 syntax, and we then remove signal part from matrix element calculation by modifying the code by hand. This is equivalent to copying matrix element computation code matrix1.f generated when we use squared order syntax:

cp template_pp_h0_tt_i/SubProcesses/P1_gg_ttx/matrix1.f template_gg_h0_tt_S_i//SubProcesses/P1_gg_ttx/.


#Running madevent
cd template_pp_h0_tt_S
./bin/madevent
launch
madspin=ON 

#Quick check of LHE content

gunzip Events/run_01/unweighted_events.lhe.gz

#You can use lhe_reader_non_decayed.c script to quickly check what is in your LHE file.

#Compilation:
#If you have TLorentz.h vector problem concerning with the root version then use root6.
find /home/muhammad/root6/ -type f -name "thisroot.sh"
source /home/muhammad/root6/build/bin/thisroot.sh
g++ -std=c++11 lhe_reader_non_decayed.c -o lhe_reader_non_decayed -I`root-config --incdir` `root-config --libs`

#g++ lhe_reader_non_decayed.c -o lhe_reader -I`root-config --incdir` `root-config --libs`

#Usage:

./lhe_reader_non_decayed unweighted_events (LHE file name without .lhe)

