# Script Created by Sayaka Miura #

from Bio import Phylo
from Bio.Phylo.Consensus import _BitString
import pandas as pd
import sys
import os

SNV = sys.argv[1]  # CloneFinder input Tree[:-19]+'.tsv'
Tree = SNV[:-4] + "snv_CloneFinder.nwk"
CloFre = Tree[:-4] + ".txt"
SV = (
    "NA"  # sys.argv[2]#Tree[:-28]+'.csv' #SV file should have the format as the example
)

OutTree = Tree[:-4] + "Tree.txt"
OutSNV = Tree[:-4] + "SNV.txt"
OutSV = Tree[:-4] + "SV.txt"
OutClone = Tree[:-4] + "CloneID.txt"

tree = Phylo.read(Tree, "newick")
CloFreCut = 0
MinVAF = 0.05

# make SNV input
SNVta = pd.read_csv(SNV, sep="\t")
print(SNVta)
Chrom = []
for i in SNVta["CHR"]:
    Chrom.append(i.replace("chr", ""))
Len = len(SNVta["CHR"])

out = {
    "#chrom": Chrom,
    "pos": SNVta["Position"],
    "desc": ["NA"] * Len,
    "normal": [0] * Len,
}
ColLs = SNVta.columns
SampLs = []
for Col in ColLs:
    if Col.find(":ref") != -1:
        SampLs.append(Col.split(":")[0])
for Samp in SampLs:
    VAF = []
    c = 0
    while c < Len:
        Tot = SNVta[Samp + ":alt"][c] + SNVta[Samp + ":ref"][c]
        if Tot == 0:
            VAF.append(0)
        else:
            VAF.append(
                1.0
                * SNVta[Samp + ":alt"][c]
                / (SNVta[Samp + ":alt"][c] + SNVta[Samp + ":ref"][c])
            )
        c += 1
    out[Samp] = VAF
out = pd.DataFrame(out)
out.to_csv(OutSNV, sep="\t", index=False)

if os.path.exists(SV) == True:
    SVta = pd.read_csv(SV, sep=",")
    print(SVta)

    ID2SV = {}
    SV2ID = {}
    ID2Samp = {}
    Rmls = []
    Len = len(SVta["Position"])
    c = 0
    while c < Len:
        if (
            SVta["#CHROM"][c].find("random") == -1
            and SVta["#CHROM"][c].find("Y") == -1
            and SVta["#CHROM"][c].find("X") == -1
            and SVta["#CHROM"][c].find("Un") == -1
        ):
            ID = SVta["ID"][c].split(":")[0]
            ALT = SVta["ALT"][c]
            REF = SVta["REF"][c]
            if ALT[-1] == REF:
                SV = (
                    ALT.replace("[", "").replace("]", "").replace(REF, "")
                    + "-"
                    + SVta["#CHROM"][c]
                    + ":"
                    + str(SVta["POS"][c])
                )
            elif ALT[0] == REF:
                SV = (
                    SVta["#CHROM"][c]
                    + ":"
                    + str(SVta["POS"][c])
                    + "-"
                    + ALT.replace("[", "").replace("]", "").replace(REF, "")
                )
            else:
                print(ALT, REF, SVta["ID"][c])
                open("a", "r").readlines()

            if SV not in SV2ID:

                SV2ID[SV] = ID
                ID2SV[ID] = SV
                ID2Samp[ID] = {}
            if ID != SV2ID[SV]:

                print(ID, SV, SV2ID[SV], ID2SV[ID])
                open("a", "r").readlines()
            for Samp in SampLs:
                Read = SVta[Samp][c].split(":")
                if float(Read[2]) < 1:
                    VAF = 0  # Rmls.append(ID)
                else:
                    VAF = 1.0 * float(Read[5]) / float(Read[2])

                    if ID2Samp[ID].get(Samp, VAF) != VAF:
                        print(
                            ID,
                            SV,
                            ID2SV[ID],
                            Samp,
                            VAF,
                            ID2Samp[ID][Samp],
                            SVta["ID"][c],
                        )
                        open("a", "r").readlines()
                ID2Samp[ID][Samp] = VAF
        c += 1
    SVnum = len(ID2Samp.keys()) + 1
    out = {"id": [], "SV": [], "Normal_GenomeCounts": []}
    for Samp in SampLs:
        out[Samp + "_GenomeCounts"] = []

    SVnum = 0
    IDin = 1
    for ID in ID2SV:
        TotVAV = 0
        Max = 0
        for Samp in SampLs:
            TotVAV += ID2Samp[ID][Samp]
            if Max < ID2Samp[ID][Samp]:
                Max = ID2Samp[ID][Samp]
        if TotVAV > 0 and MinVAF < Max and Max < 0.6:
            SVnum += 1
            out["SV"].append(ID2SV[ID])
            out["id"].append(IDin)
            IDin += 1
            for Samp in SampLs:
                out[Samp + "_GenomeCounts"].append(ID2Samp[ID][Samp])
            out["Normal_GenomeCounts"].append(0)
        else:
            print("bad SV", ID)

    print("SV count", IDin)
    out = pd.DataFrame(out)
    out.to_csv(OutSV, sep="\t", index=False)


def _bitstrs(tree):
    bitstrs = set()
    term_names = [term.name for term in tree.get_terminals()]
    term_names.sort()
    All = tree.get_terminals() + tree.get_nonterminals()
    TipID2Bit = {}
    Bit2TipID = {}
    Clone2IntBit = {}
    Dec2AncBit = {}
    for clade in All:  # tree.get_terminals():#.get_nonterminals():
        clade_term_names = [term.name for term in clade.get_terminals()]
        boolvals = [name in clade_term_names for name in term_names]
        bitstr = "".join(map(str, map(int, boolvals)))
        bitstrs.add(bitstr)
        # print (clade.get_terminals(), bitstr)
        if len(clade.get_terminals()) == 1:
            TipID2Bit[clade.get_terminals()[0].name] = bitstr
            IntLs = tree.trace(clade.get_terminals()[0].name, "Normal")
            IntBitLs = [bitstr]
            for Int in IntLs[:-2]:
                clade_term_names = [term.name for term in Int.get_terminals()]
                boolvals = [name in clade_term_names for name in term_names]
                Intbitstr = "".join(map(str, map(int, boolvals)))
                Dec2AncBit[IntBitLs[-1]] = Intbitstr
                IntBitLs.append(Intbitstr)

            Clone2IntBit[clade.get_terminals()[0].name] = IntBitLs

        Bit2TipID[bitstr] = clade.get_terminals()
    return TipID2Bit, Bit2TipID, Clone2IntBit, Dec2AncBit


TipID2Bit, Bit2TipID, Clone2IntBit, Dec2AncBit = _bitstrs(tree)
CF = pd.read_csv(CloFre, sep="\t")

CF = CF.set_index("Tumor").T
CloLs = list(CF.index)
CloC = len(CloLs)
Samp2bitVAF = {}

for Samp in SampLs:

    Sub = pd.DataFrame({"Freq": CF[Samp]})
    Sub = Sub.loc[(Sub > CloFreCut).any(axis=1)]

    CloLs = list(Sub.index)
    Len = len(CloLs)
    Bit2VAF = {}
    c = 0
    while c < Len:
        Clo = CloLs[c]
        Fre = Sub["Freq"][c] / 2

        BitLs = Clone2IntBit[Clo]
        for B in BitLs:
            Bit2VAF[B] = Bit2VAF.get(B, 0) + Fre
        c += 1
    Samp2bitVAF[Samp] = Bit2VAF

BitLs = list(Bit2TipID)

NorBit = TipID2Bit["Normal"]

MRCAbit = str(int("1" * (CloC + 1)) - int(NorBit))

BitLs.remove(NorBit)
BitLs.remove("1" * (CloC + 1))

Len = len(BitLs)

out = ["Nodes:\n"]

c = 1
Bit2ID = {}
Bit2Sbit = {}
Nout = ["NodeID\tBit\tClade\n"]
while c <= Len:
    VAFls = []
    Bit = BitLs[c - 1]
    Bit2ID[Bit] = c
    TipLs = []
    for i in Bit2TipID[Bit]:
        TipLs.append(i.name)
    Nout.append("\t".join([str(c), str(Bit), ";".join(TipLs)]) + "\n")
    Sbit = "0"
    for S in SampLs:
        VAF = Samp2bitVAF[S].get(Bit, 0)
        if VAF > 0:
            VAFls.append(VAF)
            Sbit += "1"
        else:
            Sbit += "0"
    Bit2Sbit[Bit] = Sbit

    out.append(str(c) + "\t" + Sbit + "\t" + "[ " + " ".join(map(str, VAFls)) + "]\n")
    c += 1

OutF = open(OutClone, "w")
OutF.write("".join(Nout))
OutF.close()

out.append("\n****Tree 0****\n")
out.append("0 -> " + str(Bit2ID[MRCAbit]) + "\n")

for DecB in Dec2AncBit:
    Anc = Dec2AncBit[DecB]

    out.append(str(Bit2ID[Anc]) + " -> " + str(Bit2ID[DecB]) + "\n")
out.append("Error score: 0.05\n\n")
out.append("Sample decomposition:\n")
out.append("\tSample lineage decomposition: normal\n")

out.append("\n")
for S in SampLs:
    out.append("\tSample lineage decomposition: " + S + "\n\n")

OutF = open(OutTree, "w")
OutF.write("".join(out))
OutF.close()
