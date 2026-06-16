import subprocess
import argparse
import sys
import re
########### The two parts are independent of each other#########


########### added: keep only non-overlap feature under the same BioCyc ID #########
def to_float(x):
    try:
        return float(x)
    except:
        return 0.0


def overlap(feature1, feature2):
    # Only check genomic overlap. The caller decides whether BioCyc ID is same.
    if feature1['contig_name'] != feature2['contig_name']:
        return False

    start1 = int(float(feature1['start']))
    end1 = int(float(feature1['end']))
    start2 = int(float(feature2['start']))
    end2 = int(float(feature2['end']))

    return max(start1, start2) <= min(end1, end2)


def keep_best_hit_same_biocyc(features):
    # Sort by highest bitscore first. If same bitscore, keep the earlier line.
    sorted_features = sorted(features, key=lambda x: (-to_float(x['score']), x['line_number']))

    kept = []
    removed = 0
    for feature in sorted_features:
        remove_this_feature = False
        for kept_feature in kept:
            # IMPORTANT: only remove overlap when the BioCyc ID is the same.
            # Different BioCyc IDs are allowed to overlap and are kept.
            if feature['reaction_ID'] == kept_feature['reaction_ID'] and overlap(feature, kept_feature):
                remove_this_feature = True
                removed += 1
                break

        if not remove_this_feature:
            kept.append(feature)

    # Write output in original input order, not bitscore order.
    kept = sorted(kept, key=lambda x: x['line_number'])
    print('features before filtering:', len(features))
    print('features removed because same BioCyc ID overlap:', removed)
    print('features kept:', len(kept))
    return kept
##################################################################################


if __name__== '__main__':
    name_replace = argparse.ArgumentParser()
    name_replace.add_argument('-b', help='bin_name,bin_5')
    name_replace.add_argument('-P', action='store_true',help='replace the contig name in PROKKA gtf file')
    name_replace.add_argument('-G', help='generate the GTF file from gapseq reactions file (e.g. bin.5_bin.42_Cluster.9_ce11-_ce31-_uniq_ce-all-Reactions.tbl)')
    name_replace.add_argument('-g', help='genome file name,bin.5_bin.42_Cluster.9_ce11-_ce31-_uniq_ce.fna')
    name_replace.add_argument('-T', help='enerate the GTF file from gapseq trans file,bin.5_bin.42_Cluster.9_ce11-_ce31-_uniq_ce-Transporter.tbl')
    name_replace.add_argument('-o', default='output.gtf',help='output file')
    args = vars(name_replace.parse_args())

    output=args['o']
############# first part: for the GTF file generated from PROKKA gff file #################
    ### since this script use linux command, it can not run in the windows (local) for fitst part
    ### input the bin name and genome file name, and the target file (gtf file with incorrect contig name that can not be used in htseq-count)
    ### output a new gtf file with correct contig name based on genome file
    ### replace_contig_name_prokka_20230815.py -b bin_5 -g bin.5_bin.42_Cluster.9_ce11-_ce31-_uniq_ce -o replaced_right_name.gtf
    ### extract the contig name from genome file and store in contig_name_list
    if '-P' in sys.argv:
        command_contig_name = "grep '>' "+args['g']+" | sed 's/>//g'"

        result = subprocess.run(command_contig_name, stdout=subprocess.PIPE, shell=True)

        contig_name_str = result.stdout.decode('utf-8') ## use \n as the separator
        contig_name_list=contig_name_str.strip('\n').split('\n')

        #### extract the prokka contig name and store in prokka_name_list
        #command_prokka="awk '/^##sequence-region/{print$2}' bin_5.gff|awk '!seen[$0]++'"
        command_prokka="awk '/^##sequence-region/{print$2}' PROKKA/"+args['b']+".gff|awk '!seen[$0]++'"
        result_prokka=subprocess.run(command_prokka, stdout=subprocess.PIPE, shell=True)

        prokka_name_str=result_prokka.stdout.decode('utf-8')
        prokka_name_list=prokka_name_str.strip('\n').split('\n')
        ### put both in the dict and read the dict simultenously and output a new file with replaced line
        replacement_map=dict(zip(prokka_name_list,contig_name_list))

        with open(args['b']+'.gtf','r') as infile, open(output,'w') as outfile:
            for line in infile:
                for old,new in replacement_map.items():
                    line = line.replace(old, new)
                outfile.write(line)

################### second part: for gapseq reactions file ###############
###### if the transcriptomic data is used to study the expression of gapseq reactions
##### generate the gtf file from -all-Reactions.tbl file #######
### gtf file in form: contig_name,data_source,feature_type,start and end position,bitscore,strand,frame,gene_attributes
    if '-G' in sys.argv:
        reactions_file=open(args['G']).readlines()
        source='gapseq'
        feature_type='CDS'
        features=[]

        for line_number, line in enumerate(reactions_file):
            if not line.startswith('#'): #remove the first line in the gapseq reactions file
                try:
                    reactions_status =line.strip('\n').split('\t')[13]### consider the good blast and bad blast reactions
                except:
                    print('a line does not have 14th element')
                    continue
                if reactions_status=='good_blast':

                    contig_name=line.strip('\n').split('\t')[9]
                    score=line.strip('\n').split('\t')[7]
                    reaction_ID=line.strip('\n').split('\t')[0]
                    rxn_name=line.strip('\n').split('\t')[1]
                    subunit=line.strip('\n').split('\t')[16]
                    attributes='gene_id '+reaction_ID+'_'+reactions_status+'_'+subunit+';Name='+rxn_name+';Note='+subunit
                ##### judge the direction of CDS #####
                    first_point=line.strip('\n').split('\t')[10]
                    last_point=line.strip('\n').split('\t')[11]
                    if float(first_point)>float(last_point):
                        strand='-'
                        start=last_point
                        end=first_point
                    else:
                        strand='+'
                        start=first_point
                        end=last_point
        ############  store the new line, then filter overlaps under same BioCyc ID #########
                    new_line='\t'.join([contig_name,source,feature_type,start,end,score,strand,'.',attributes])+'\n'
                    features.append({'contig_name':contig_name,
                                     'start':start,
                                     'end':end,
                                     'score':score,
                                     'reaction_ID':reaction_ID,
                                     'line_number':line_number,
                                     'new_line':new_line})

        features=keep_best_hit_same_biocyc(features)
        with open(output, 'w') as outfile:
            for feature in features:
                outfile.write(feature['new_line'])


##### the third part is for gapseq transporter file, most of this is similar with last part.
    if '-T' in sys.argv:
        reactions_file=open(args['T']).readlines()
        source='gapseq_transporter'
        feature_type='CDS'
        features=[]

        for line_number, line in enumerate(reactions_file[3:]): # remove the first line in the gapseq reactions file

            contig_name = line.strip('\n').split('\t')[10]
            score = line.strip('\n').split('\t')[8]
            reaction_ID = line.strip('\n').split('\t')[1]
            rxn_name = line.strip('\n').split('\t')[2]

            compound_ID_EX = line.strip('\n').split('\t')[3]
            compound_ID_match=re.search(r'_(.*?)_',compound_ID_EX)
            compound_ID=compound_ID_match.group(1)

            attributes = 'gene_id ' + reaction_ID+'_'+rxn_name + '_'+compound_ID+';Name=' + rxn_name + ';Note=' + compound_ID

            ##### judge the direction of CDS #####
            first_point=line.strip('\n').split('\t')[11]
            last_point=line.strip('\n').split('\t')[12]
            if float(first_point)>float(last_point):
                strand='-'
                start=last_point
                end=first_point
            else:
                strand='+'
                start=first_point
                end=last_point
    ############  store the new line, then filter overlaps under same BioCyc ID #########
            new_line='\t'.join([contig_name,source,feature_type,start,end,score,strand,'.',attributes])+'\n'
            features.append({'contig_name':contig_name,
                             'start':start,
                             'end':end,
                             'score':score,
                             'reaction_ID':reaction_ID,
                             'line_number':line_number,
                             'new_line':new_line})

        features=keep_best_hit_same_biocyc(features)
        with open(output, 'w') as outfile:
            for feature in features:
                outfile.write(feature['new_line'])
