import os
import datetime
import pandas as pd
from outbreak_data import outbreak_data

cryptics_dir = 'cryptic'
max_gisaid_sequences = 30

def sort_func(x):
    if 'DEL' in x:
        if '/' in x:
            return int(x.split('DEL')[1].split('/')[0])
        else:
            return int(x.split('DEL')[1])
    else:
        return int(x.split(':')[1][1:-1])

def process_cryptic(mutations):

    mutations = [mut.replace("'", "").strip() for mut in mutations]
    mutations = sorted(mutations, key=sort_func)

    return mutations

common_muts = []
rare_muts = []
def find_related_lineages(covariants):
    covariants = covariants.split(',')
    # Find the rare mutation
    for mut in covariants:
        if mut in common_muts or mut in rare_muts:
            continue

        # TODO: add location to query
        query = f"mutations={mut}"
        try:
            results = outbreak_data.get_outbreak_data(
                "genomics/mutations-by-lineage", argstring=query
            )["results"]
            common_muts.append(mut)
        except NameError as e:  # No clinical results found
            rare_muts.append(mut)
            continue

    # remove rare muts and find the corresponding lineage
    ancestral_cluster = [mut for mut in covariants if mut not in rare_muts]
    
    mutations = ','.join(ancestral_cluster)
    query = f"mutations={mutations}"
    try:
        results = outbreak_data.get_outbreak_data(
            "genomics/mutations-by-lineage", argstring=query
        )["results"]
        results = list(results[mutations])
    except NameError as e:  # No clinical results found
        print(query)
        print('rare muts', rare_muts)
        return []
    lineages = [dict['pangolin_lineage'] for dict in results]
    lineages = list(set(lineages))
    return f'[{",".join(lineages)}]'

cryptics_list = pd.DataFrame(columns=['Mutation Cluster', '# Detections', 'First Detection', 'Last Detection', 'Samples Detected'])

for file in os.listdir(cryptics_dir):
    df = pd.read_csv(os.path.join(cryptics_dir, file), sep='\t')
    
    df['Covariants'] = df['Covariants'].apply(lambda x: x.strip('][').split(','))
    df['Covariants'] = df['Covariants'].apply(process_cryptic)
    
    # Drop rows with no covariants
    df = df.dropna(subset=['Covariants'])
    
    if not df.empty:
        cryptic_muts = df['Covariants'].tolist()

        for mut in cryptic_muts:
            date = datetime.datetime.strptime(file.split('_')[0], '%Y-%m-%d').date()
            sewershed = file.split('_')[1].split('.')[0]
            mut = ','.join(mut)
            if mut not in cryptics_list['Mutation Cluster'].tolist():
                cryptics_list.loc[len(cryptics_list.index)] = [mut, 1, date, date, f'{sewershed}({date})']
            else:
                cryptics_list.loc[cryptics_list['Mutation Cluster'] == mut,'# Detections'] += 1
                
                if date < cryptics_list.loc[cryptics_list['Mutation Cluster'] == mut,'First Detection'].item():
                    cryptics_list.loc[cryptics_list['Mutation Cluster'] == mut,'First Detection'] = date
                if date > cryptics_list.loc[cryptics_list['Mutation Cluster'] == mut,'Last Detection'].item():
                    cryptics_list.loc[cryptics_list['Mutation Cluster'] == mut,'Last Detection'] = date

                cryptics_list.loc[cryptics_list['Mutation Cluster'] == mut,'Samples Detected'] += f' {sewershed}({date})'

#cryptics_list['Related Lineage(s)'] = cryptics_list['Mutation Cluster'].apply(find_related_lineages)
cryptics_list['Samples Detected'] = cryptics_list['Samples Detected'].apply(lambda x: x.split(' '))
cryptics_list['Samples Detected'] = cryptics_list['Samples Detected'].apply(sorted, key=lambda x: datetime.datetime.strptime(x.split('(')[1].split(')')[0], '%Y-%m-%d').date())

cryptics_list = cryptics_list[cryptics_list['# Detections'] > 1]
cryptics_list = cryptics_list.sort_values(by=['# Detections'], ascending=False)

cryptics_list.to_csv('cryptics_list.tsv', sep='\t', index=False)

