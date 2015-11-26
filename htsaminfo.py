import pandas as pd

def parse_opt(l_options):
    gene = "__"
    for opt in l_options:
        f_opt = opt.split(":")
        if f_opt[0] == "XF" and f_opt[1] == "Z":
            gene = f_opt[2]
            break

    return gene


class HTSeqSamInfo(object):
    def __init__(self, input_sam, geneList, erccList):
        self.insam = input_sam
        self.genList = {}
        self.genList['gene'] = self.__read_list(geneList)
        self.genList['ercc'] = self.__read_list(erccList)
        self.genUMIcnt = {"gene":{}, "ercc":{} }
        self.UMI_info = {"gene": [], "UMI": [], "type":[]}

    def s1_SamToUMI(self):
        with open(self.insam, "r") as f_sam:
            for line in f_sam:
                if line[0] == "@":
                    continue
                line = line.strip()
                f = line.split()
                UMI = f[0][0:8]
                l_options = f[11:]
                gene = parse_opt(l_options)
                if gene[0:2] == "__" or gene[0:3].upper() == "MIR" or  gene[0:5].upper() == "SNORD":
                    continue

                tag = "gene"
                if gene[0:5].upper() == "ERCC-" or gene[0:4].upper() == "RGC-":
                    tag = "ercc"

                if gene not in self.genUMIcnt[tag]:
                    self.genUMIcnt[tag][gene] = {}

                if UMI  not in self.genUMIcnt[tag][gene]:
                    self.genUMIcnt[tag][gene][UMI] = 0

                self.genUMIcnt[tag][gene][UMI] += 1


    def s2_outUMIList(self):
        for tag in ['gene', 'ercc']:
            for gene in self.genList[tag]:
                UMI = 0
                if gene in self.genUMIcnt[tag]:
                    UMI = len(self.genUMIcnt[tag][gene])
                self.UMI_info['gene'].append(gene)
                self.UMI_info['UMI'].append(UMI)
                self.UMI_info['type'].append(tag)

        self.df_UMIInfo = pd.DataFrame(self.UMI_info)
        outFile = "%s.UMI_gene.xls" % (".".join(self.insam.split(".")[:-1]))
        self.df_UMIInfo.to_csv(outFile, sep="\t", index=None)
        
    def s3_UMI2TPM(self):
        self.df_UMIInfo['TPM'] = self.df_UMIInfo['UMI'] * 1*10**6 / \
                                 self.df_UMIInfo['UMI'].sum()
        
        outFile = "%s.UMI_withTPM.xls" % (".".join(self.insam.split(".")[:-1]))
        self.df_UMIInfo.to_csv(outFile, sep="\t", index=None)

    
    """
        Private Functions:
    """
    def __read_list(self, infile):
        l_list = []
        with open(infile, "r") as f_infile:
            for line in f_infile:
                line = line.strip()
                l_list.append(line)
        return l_list
                