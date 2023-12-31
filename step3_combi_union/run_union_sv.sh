#!/usr/bin/sh
#SBATCH -A MST109178        # Account name/project number
#SBATCH -J union_HG002    # Job name
#SBATCH -p ngs186G           # Partition Name 等同PBS裡面的 -q Queue name
#SBATCH -c 28              # 使用的數 請參考Queue資源設定
#SBATCH --mem=186G           # 使用的記憶體量 請參考Queue資源設定
#SBATCH -o out.union_HG002.log          # Path to the standard output file
#SBATCH -e err.union_HG002.log          # Path to the standard error ouput file
#SBATCH --mail-user=none    # email
#SBATCH --mail-type=BEGIN,END              # 指定送出email時機 可為NONE, BEGIN, END, FAIL, REQUEUE, AL

module load biology/Python/3.11.0

#inputdir= sys.argv[1]
#outputdir= sys.argv[2]
#mergefile= sys.argv[3]
#firstfile= sys.argv[4]
#allfile= sys.argv[5]
#statfile= sys.argv[6]

python ./all_in_one_sv0.py  HG002.mdg union.HG002.mdg mergedsv.mdg union.first.mdg union.all.mdg union.stat.mdg
python ./all_in_one_sv0.py  HG002.mdgls union.HG002.mdgls mergedsv.mdgls union.first.mdgls union.all.mdgls union.stat.mdgls
