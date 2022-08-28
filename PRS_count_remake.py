import pandas as pd
import io, sys, os
from numpy import log as ln
def read_vcf(path):
    with open(path, 'r') as f:
        lines = [l for l in f if not l.startswith('##')]
    return pd.read_csv(
        io.StringIO(''.join(lines)),
        dtype={'#CHROM': str, 'POS': int, 'ID': str, 'REF': str, 'ALT': str,
               'QUAL': str, 'FILTER': str, 'INFO': str},
        sep='\t'
    ).rename(columns={'#CHROM': 'CHROM'})

paper_HBOC_data = pd.read_csv("41467_2020_17680_MOESM4_ESM_transfered_to_38.csv")
paper_HBOC_data.info()
position = paper_HBOC_data["Position"].to_list()
Chr = paper_HBOC_data["Chr."].to_list()
Effect_size = paper_HBOC_data["Chinese OR"].to_list()
allele = paper_HBOC_data["Allele"].to_list()

data_dir = sys.argv[1]
file_name = sys.argv[2]
# print(file_name)
vcf = read_vcf(data_dir)
vcf_CHROM = vcf["CHROM"].to_list()
vcf_Position = vcf["POS"].to_list()
csv_rows = []
final = []
score = []
for i in range(0,len(Chr)):
    temp = []
    # print(position[i])
    # print(Chr[i])
    count = len(final)
    for j in range(0,len(vcf_CHROM)):
        temp_allele = vcf["REF"][j] + "/" +vcf["ALT"][j]
        if str(vcf_CHROM[j].replace("chr","")) == str(Chr[i]) and str(vcf_Position[j]) == str(position[i]) and str(temp_allele) == str(allele[i]):
            # print(str(vcf_CHROM[j].replace("chr","")))
            # print(str(Chr[i]))
            # print(vcf[file_name][i].split(":")[0])
            if vcf[file_name][i].split(":")[0] == "0/1":
                dosage = 1
            elif vcf[file_name][i].split(":")[0] == "1/1":
                dosage = 2
            elif vcf[file_name][i].split(":")[0] == "0/0":
                dosage = 0
            # print(Effect_size[i])
            # print(ln(Effect_size[i]))
            temp_fianl = ln(Effect_size[i]) * dosage
            final.append(temp_fianl)
            print(f"=====>> {i} <<======")
            temp.append([paper_HBOC_data["SNP"][i], vcf_CHROM[j], position[i], paper_HBOC_data["Allele"][i], dosage, Effect_size[i], ln(Effect_size[i]), temp_fianl])
            # print(temp)
            score.append(float(temp_fianl))
            csv_rows.append(temp[0])
    new_count  = len(final)
    if count == new_count:
        temp.append([paper_HBOC_data["SNP"][i], f"chr{Chr[i]}", position[i], paper_HBOC_data["Allele"][i], "null", Effect_size[i], ln(Effect_size[i]), "null"])
        # print(temp)
        csv_rows.append(temp[0])
    cols = ["SNP ID", "CHR", "BP", "Allele", "chinese OR", "Effect Size", "sample 1", "score"]
    import csv
    with open(f"/home/Digi118078/0808/PRS_score/{file_name}_result.csv","w") as f:
        write = csv.writer(f)
        write.writerow(cols)
        write.writerows(csv_rows)
with open("/home/Digi118078/0808/PRS_score/score.txt","w") as score_file:
    score_file.write(f"{file_name}\t{sum(count)}\t{len(final)}\n")
# print("+++++++++++++++")
# print(sum(count))
# print(len(final))





# cols = ["SNP ID", "CHR", "BP", "Allele", "chinese OR", "Effect Size", "sample 1", "score"]
# import csv
# # rows = [f_SnpID, f_Chr, f_Bp, f_Allele, f_sample_1, f_chi_OR, f_EffectSize, f_score]
# with open(f"result_{name}.csv","w") as f:
#     write = csv.writer(f)
#     write.writerow(cols)
#     write.writerows(csv_rows)

# # for i in range(0,len(vcf_CHROM)):
# #     print(vcf_CHROM[i].replace("chr",""))  #1
# #     print(vcf["file_name"][i].split(":")[0])   #0/1
# #     if str(vcf_CHROM[i].replace("chr","")) == '1':

# #     break
# #     #     print(vcf["file_name"][i])
# #     #     break

# '''
# SNP ID | CHR | BP | Allele | chinese OR | Effect Size | sample 1 | score

# f_SnpID = []
# f_Chr = []
# f_Bp = []
# f_Allele = []
# f_sample_1 = []
# f_chi_OR = []
# f_EffectSize = []
# f_score = []
# '''