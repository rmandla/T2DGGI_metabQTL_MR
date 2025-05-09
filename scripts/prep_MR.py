import pandas as pd
import sys, subprocess, argparse, os

class customError(Exception):
    def __init__(self, message):
        super().__init__(message)

def run_clumping(sst,ref_path,exposure,chrom_col,bp_col,P_col,ea_col,nea_col,isLogP=False,output_header='',dataset=None,plink='plink',snps=None,pthresh=5e-8,rsid_mappings=None):
    sst_df = pd.read_table(sst)
    if exposure.lower() == 'mqtl':
        if output_header != '':
            output_header += '.'
        output = f'{output_header}{exposure}.{dataset}'
        if chrom_col == '':
            raise customError('chromosome column name needs to be specified')
        if bp_col == '':
            raise customError('base pair column name needs to be specified')
        if ea_col == '':
            raise customError('effect allele column name needs to be specified')
        if nea_col == '':
            raise customError('non-effect allele column name needs to be specified')
        if P_col == '':
            raise customError('p-value column name needs to be specified')
    elif exposure.lower() == 'gwas':
        if output_header != '':
            output_header += '.'
        output = f'{output_header}{exposure}.T2DGGI'
        if chrom_col == '':
            chrom_col='CHR'
        if bp_col == '':
            bp_col='BP'
        if ea_col == '':
            ea_col='Allele1'
        if nea_col == '':
            nea_col='Allele2'
        if P_col == '':
            P_col='P'
    else:
        raise(f'ERROR: {exposure} not recognized')

    if type(snps) == str:
        snps_df = pd.read_table(snps,header=None)
        if 'RSID' in sst_df.columns:
            sst_df = sst_df[sst_df['RSID'].isin(snps_df[0])]
        #elif 'rsids' in sst_df.columns:
        #    sst_df = sst_df[sst_df['rsids'].isin(snps_df[0])]
        elif type(rsid_mappings) == str:
            rsid_mappings_df = pd.read_table(rsid_mappings)
            sst_df['SNP_ID1'] = sst_df[chrom_col].astype(str)+':'+sst_df[bp_col].astype(str)+':'+sst_df[ea_col]+':'+sst_df[nea_col]
            sst_df['SNP_ID2'] = sst_df[chrom_col].astype(str)+':'+sst_df[bp_col].astype(str)+':'+sst_df[nea_col]+':'+sst_df[ea_col]
            temp_sst_df1 = sst_df.merge(rsid_mappings_df,left_on='SNP_ID1',right_on='SNP').drop(columns=['SNP_ID1'])
            temp_sst_df2 = sst_df.merge(rsid_mappings_df,left_on='SNP_ID2',right_on='SNP').drop(columns=['SNP_ID2'])
            sst_df = pd.concat([temp_sst_df1,temp_sst_df2])
            sst_df = sst_df[sst_df['RSID'].isin(snps_df[0])]
        else:
            raise(f'ERROR: RSID column not found in input file and no alternative RSID mappings were provided')

    sst_df['CHR'] = sst_df[chrom_col]
    sst_df = sst_df[sst_df['CHR'].isin([str(i) for i in range(1,23)]+[i for i in range(1,23)])]
    sst_df['CHR'] = sst_df['CHR'].astype(int)
    sst_df['BP'] = sst_df[bp_col]
    if isLogP:
        sst_df['P'] = 10**-sst_df[P_col]
    else:
        sst_df['P'] = sst_df[P_col]
    sst_df['SNP1'] = sst_df['CHR'].astype(str)+':'+sst_df['BP'].astype(str)+':'+sst_df[ea_col].str.upper()+':'+sst_df[nea_col].str.upper()
    sst_df['SNP2'] = sst_df['CHR'].astype(str)+':'+sst_df['BP'].astype(str)+':'+sst_df[nea_col].str.upper()+':'+sst_df[ea_col].str.upper()
    sst_df1 = sst_df.copy()
    sst_df1['SNP'] = sst_df1['SNP1']
    sst_df2 = sst_df.copy()
    sst_df2['SNP'] = sst_df2['SNP2']
    sst_df = pd.concat([sst_df1,sst_df2])
    if '/' in output:
        output_dir = '/'.join(output.split('/')[:-1])+'/'
        output_file = output.split('/')[-1]
    else:
        output_dir = './'
        output_file = output
    if output+".ALL.clumped" in os.listdir(output_dir):
        subprocess.run(f'rm {output}.ALL.clumped',shell=True,check=True)
    output_files = []
    for chrom in range(1,23):
        tsst_df = sst_df[sst_df['CHR']==chrom]
        tsst_df[['CHR','SNP','BP','P']].to_csv('temp_forclump.txt',sep='\t',index=None)
        ref_file = f'{ref_path}/EUR_1KG_chr{chrom}'
        subprocess.run(f'{plink} --bfile {ref_file} --clump temp_forclump.txt --clump-p1 {pthresh} --clump-kb 10000 --clump-r2 0.001 --out {output}.chr{chrom}',shell=True,check=True)

        if f'{output_file}.chr{chrom}.clumped' in os.listdir(output_dir):
            output_files.append(f'{output_dir}{output_file}.chr{chrom}.clumped')
    if len(output_files)>0:
        for of in output_files:
            tdf = pd.read_table(of,delim_whitespace=True)[['CHR','SNP','BP']]
            if of == output_files[0]:
                ofs = tdf.copy()
            else:
                ofs = pd.concat([ofs,tdf])
        ofs.to_csv(f'{output_dir}{output_file}.ALL.clumped',sep='\t',index=None)

def prep_GWAS_data(gwas_path,metabname,dataset,output_header,clumped_snps=None):
    metal = pd.read_table(gwas_path)

    metal = metal[['CHR','BP','RSID','Allele1','Allele2','P','Effect','StdErr','Freq1','N','MarkerName']]
    metal.columns = ['chr','pos','rsid','eAllele','oAllele','P','Effect','SE','AF','N','snpid']
    metal['eAllele'] = metal['eAllele'].str.upper()
    metal['oAllele'] = metal['oAllele'].str.upper()

    if type(clumped_snps) == str:
        snps = pd.read_table(clumped_snps,delim_whitespace=True)
        snps['snpid'] = 'chr'+snps['CHR'].astype(str)+':'+snps['BP'].astype(str)
        metal = metal[metal['snpid'].isin(snps['snpid'])]
    metal.drop_duplicates().sort_values('rsid').drop(columns=['snpid']).to_csv(f'{output_header}_{metabname}_{dataset}_GWAS_harmonized.txt',sep='\t',index=None)
    N = metal.sort_values('N',ascending=False)['N'].to_list()[0]
    return(N)

def prep_mQTL_data(metabname,dataset,mQTL_path,output_header,chrom_col,bp_col,beta_col,se_col,af_col,P_col,N_col,ea_col,nea_col,isLogP=False,clumped_snps=None):
    gwas = pd.read_table(f'{output_header}_{metabname}_{dataset}_GWAS_harmonized.txt')
    mqtl = pd.read_table(mQTL_path)

    gwas['SNP'] = gwas['chr'].astype(str)+':'+gwas['pos'].astype(str)+':'+gwas['eAllele']+':'+gwas['oAllele']
    rsid = gwas[['rsid','SNP']]
    mqtl['Chrom'] = mqtl[chrom_col].astype(str).str.split('chr',expand=True)[1]
    mqtl['SNP1'] = mqtl['Chrom']+':'+mqtl[bp_col].astype(str)+':'+mqtl[ea_col]+':'+mqtl[nea_col]
    mqtl['SNP2'] = mqtl['Chrom']+':'+mqtl[bp_col].astype(str)+':'+mqtl[nea_col]+':'+mqtl[ea_col]

    mqtl1 = mqtl[mqtl['SNP1'].isin(gwas['SNP'])]
    mqtl1['SNP'] = mqtl1['SNP1']
    mqtl2 = mqtl[mqtl['SNP2'].isin(gwas['SNP'])]
    mqtl2['SNP'] = mqtl2['SNP2']

    mqtl1 = mqtl1.merge(rsid,left_on='SNP1',right_on='SNP')
    mqtl2 = mqtl2.merge(rsid,left_on='SNP2',right_on='SNP')

    mqtl = pd.concat([mqtl1,mqtl2])
    if N_col not in mqtl.columns:
        try:
            N = float(N_col)
            mqtl['N'] = N
            N_col = 'N'
        except:
            raise customError(f'{N_col} not found in mQTL dataset')
    mqtl = mqtl[['Chrom',bp_col,'rsid',ea_col,nea_col,P_col,beta_col,se_col,af_col,N_col]]
    if isLogP:
        mqtl[P_col] = 10**(-mqtl[P_col])
    if np.max(mqtl[af_col]>0.5):
        mqtl[af_col][mqtl[af_col]>0.5] = 1-mqtl[af_col]

    mqtl.columns = ['chr','pos','rsid','eAllele','oAllele','P','Effect','SE','AF','N']
    if type(clumped_snps) == str:
        snps = pd.read_table(clumped_snps,delim_whitespace=True)
        snps['snpid'] = snps['CHR'].astype(str)+':'+snps['BP'].astype(str)
        mqtl['snpid'] = mqtl['chr'].astype(str)+':'+mqtl['pos'].astype(str)
        mqtl = mqtl[mqtl['snpid'].isin(snps['snpid'])]
        mqtl = mqtl.drop(columns=['snpid'])
    mqtl.sort_values('rsid').to_csv(f'{output_header}_{metabname}_{dataset}_mQTL_harmonized.txt',sep='\t',index=None)
    N = mqtl.sort_values('N',ascending=False)['N'].to_list()[0]
    return(N)

def main():
    parser = argparse.ArgumentParser()
    #subparser = parser.add_subparsers()

    subparsers = parser.add_subparsers(title='Commands', dest='command')
    # use dispatch pattern to invoke method with same name

    #parser.add_argument('pull_ld')
    clump = subparsers.add_parser('clump')
    clump.add_argument("-s", "--sst",dest='sst')
    clump.add_argument('-r','--ref_path',dest='ref_path')
    clump.add_argument('-o','--output_header',dest='output_header',default='')
    clump.add_argument('-e','--exposure',dest='exposure')
    clump.add_argument('-d','--dataset',dest='dataset',default=None)
    clump.add_argument('-p','--plink_path',dest='plink_path',default='plink')
    clump.add_argument('-pt','--p-thresh',dest='pthresh',default=5e-8)
    clump.add_argument('-snps','--snps',dest='snps',default=None)
    clump.add_argument('-rsids','--rsids',dest='rsid_mappings',default=None)
    clump.add_argument('-c','--chrom_col',dest='chrom_col',default='')
    clump.add_argument('-b','--bp_col',dest='bp_col',default='')
    clump.add_argument('-pv','--P_col',dest='P_col',default='')
    clump.add_argument('-ea','--ea_col',dest='ea_col',default='')
    clump.add_argument('-nea','--nea_col',dest='nea_col',default='')
    clump.add_argument('-isLogP','--isLogP',dest='isLogP',default=False)

    #parser.add_argument('-c','--compare',dest='compare')
    mr = subparsers.add_parser('run_mr')
    mr.add_argument('-s','--sst',dest='gwas_path')
    mr.add_argument("-p", "--protein_name",dest='metabname')
    mr.add_argument('-d', '--dataset',dest='dataset')
    mr.add_argument('-ex','--exposure',dest='exposure')
    mr.add_argument('-oc','--outcome',dest='outcome')
    mr.add_argument('-dpath','--dataset_path',dest='mQTL_path')
    mr.add_argument('-cs','--clumped_snps',dest='clumped_snps',default=None)
    mr.add_argument("-o", "--output",dest = "output_header")
    mr.add_argument('-c','--chrom_col',dest='chrom_col',default=None)
    mr.add_argument('-b','--bp_col',dest='bp_col',default=None)
    mr.add_argument('-be','-beta_col',dest='beta_col',default=None)
    mr.add_argument('-se','--se_col',dest='se_col',default=False)
    mr.add_argument('-pv','--P_col',dest='P_col',default=None)
    mr.add_argument('-ea','--ea_col',dest='ea_col',default=None)
    mr.add_argument('-nea','--nea_col',dest='nea_col',default=None)
    mr.add_argument('-isLogP','--isLogP',dest='isLogP',default=False)
    mr.add_argument('-N','--N_col',dest='N_col',default=False)
    mr.add_argument('-a','--af_col',dest='af_col',default=False)

    args = parser.parse_args()
    if args.isLogP != False:
        args.isLogP = True

    if args.command == 'clump':
        print('Clumping variants')
        run_clumping(sst=args.sst,ref_path=args.ref_path,exposure=args.exposure,chrom_col=args.chrom_col,bp_col=args.bp_col,P_col=args.P_col,ea_col=args.ea_col,nea_col=args.nea_col,output_header=args.output_header,dataset=args.dataset,plink=args.plink_path,snps=args.snps,pthresh=float(args.pthresh),rsid_mappings=args.rsid_mappings)
    elif args.command == 'run_mr':
        print('Running MR')
        gwas_n = prep_GWAS_data(gwas_path=args.gwas_path,metabname=args.metabname,dataset=args.dataset,output_header=args.output_header,clumped_snps=args.clumped_snps)
        mqtl_n = prep_mQTL_data(metabname=args.metabname,N_col=args.N_col,af_col=args.af_col,chrom_col=args.chrom_col,bp_col=args.bp_col,beta_col=args.beta_col,se_col=args.se_col,isLogP=args.isLogP,P_col=args.P_col,ea_col=args.ea_col,nea_col=args.nea_col,dataset=args.dataset,mQTL_path=args.mQTL_path,output_header=args.output_header,clumped_snps=args.clumped_snps)
        exp = args.exposure.upper()
        oc = args.outcome.upper()
        if exp == 'GWAS' and oc == 'MQTL':
            oc = 'mQTL'
            N_outcome = mqtl_n
            N_exposure = gwas_n
        elif exp == 'MQTL' and oc == 'GWAS':
            exp = 'mQTL'
            N_outcome = gwas_n
            N_exposure = mqtl_n
        else:
            raise('exposure and outcome must either be GWAS or MQTL')
        file_path = os.path.realpath(__file__)
        file_path = '/'.join(file_path.split('/')[:-1])
        for filter in ['fstat','steig','both']:
            subprocess.run(f'Rscript {file_path}/run_MR.R {args.metabname} {args.dataset} {exp} {oc} {N_exposure} {N_outcome} {args.output_header} {filter}',shell=True)
    else:
        # If an unknown command is provided, show the help message
        parser.print_help()

if __name__ == "__main__":
    main()
