process pystats {
    input:
        val mypath
    output:
        stdout
        //val mypath
        //path "pyoutputs.txt", emit: pyoutputs
        
    $/
    #!/usr/bin/env python3
    import subprocess
    
    items = "${mypath}".strip().split("/")
    #print(items[-1])
    filepath1 = "${mypath}"+"/alignment/"+items[-1]+".coverage.txt"
    #print(filepath1)
    with open(filepath1, 'r') as cov_report:
        header = cov_report.readline()
        header = header.rstrip()
        stats = cov_report.readline()
        stats = stats.rstrip()
        stats = stats.split()
        ref_name = stats[0]
        #print(ref_name)
        start = stats[1]
        end = stats[2]
        reads_mapped = stats[3]
        cov_bases = stats[4]
        cov = stats[5]
        depth = stats[6]
        baseq = stats[7]
        #print(reads_mapped)
        mapq = stats[8]
        
    #Get number of raw reads
    proc_1 = subprocess.run('zcat ' + "${mypath}/" + items[-1] + '_1_humanclean.fastq.gz | wc -l', shell=True, capture_output=True, text=True, check=True)
    wc_out_1 = proc_1.stdout.rstrip()
    reads_1 = int(wc_out_1) / 4
    proc_2 = subprocess.run('zcat ' + "${mypath}/" + items[-1] + '_2_humanclean.fastq.gz | wc -l', shell=True, capture_output=True, text=True, check=True)
    wc_out_2 = proc_2.stdout.rstrip()
    reads_2 = int(wc_out_2) / 4
    raw_reads = reads_1 + reads_2
    raw_reads = int(raw_reads)

    #Get number of clean reads
    proc_c1x = subprocess.run('zcat ' + "${mypath}/" + items[-1] + '_1.fq.gz | wc -l', shell=True, capture_output=True, text=True, check=True)
    wc_out_c1x = proc_c1x.stdout.rstrip()
    reads_c1x = int(wc_out_c1x) / 4
    proc_c2x = subprocess.run('zcat ' + "${mypath}/" + items[-1] + '_2.fq.gz | wc -l', shell=True, capture_output=True, text=True, check=True)
    wc_out_c2x = proc_c2x.stdout.rstrip()
    reads_c2x = int(wc_out_c2x) / 4
    clean_reads = reads_c1x + reads_c2x
    clean_reads = int(clean_reads)
    #print(clean_reads)
    
    #Get percentage of mapped reads/clean reads
    percent_map = "%0.4f"%((int(reads_mapped)/int(clean_reads))*100)
    #print(percent_map)
    
    #Gather QC metrics for consensus assembly
    filepath2 = "${mypath}"+"/assembly/"+items[-1]+".consensus.fa"
    with open(filepath2, 'r') as assem:
        header = assem.readline()
        header = header.rstrip()
        bases = assem.readline()
        bases = bases.rstrip()
        num_bases = len(bases)
        ns = bases.count('N')
        called = num_bases - ns
        pg = "%0.4f"%((called/int(end))*100)
        #print(called)
        #print(end)
        #print(pg)
    #Rename header in fasta to just sample name
    subprocess.run("sed -i \'s/^>.*/>"+items[-1]+"/\' "+filepath2, shell=True, check=True)
    #print("sed -i \'s/^>.*/>"+items[-1]+"/\' "+filepath2)
    
    #QC flag
    pg_flag = ''
    dp_flag = ''
    qc_flag = ''
    if float(pg) < 79.5:
        pg_flag = 'FAIL: Percent genome < 80%'
        qc_flag = qc_flag + pg_flag
    else:
        if float(depth) < 100:
            dp_flag = 'FAIL: Mean read depth < 100x'
            qc_flag = qc_flag + dp_flag
        if qc_flag == '':
            qc_flag = qc_flag + 'PASS'
    #print(qc_flag)
    
    if qc_flag == 'PASS':
        subprocess.run("cp "+filepath2+" "+"${params.output}"+"/assemblies/", shell=True, check=True)   
        subprocess.run('cp ' + "${mypath}" + '/variants/' + items[-1] + '.variants.tsv ' + "${params.output}"+'/variants/', shell=True, check=True)
    
        #Run VADR
        out_log = open("${mypath}/"+items[-1]+'.out', 'w')
        err_log = open("${mypath}/"+items[-1]+'.err', 'w')
        #subprocess.run("singularity exec -B "+"${mypath}"+"/assembly"+":/data /apps/staphb-toolkit/containers/vadr_1.3.sif /opt/vadr/vadr/miniscripts/fasta-trim-terminal-ambigs.pl --minlen 50 --maxlen 30000 /data/" + items[-1]+".consensus.fa > " + "${mypath}"+"/assembly/"+items[-1]+".trimmed.fasta", shell=True, stdout=out_log, stderr=err_log, check=True)
        subprocess.run("singularity exec -B "+"${mypath}"+"/assembly"+":/data docker://staphb/vadr:1.3 /opt/vadr/vadr/miniscripts/fasta-trim-terminal-ambigs.pl --minlen 50 --maxlen 30000 /data/" + items[-1]+".consensus.fa > " + "${mypath}"+"/assembly/"+items[-1]+".trimmed.fasta", shell=True, stdout=out_log, stderr=err_log, check=True)

        #subprocess.run("singularity exec -B "+"${mypath}"+"/assembly:/data /apps/staphb-toolkit/containers/vadr_1.3.sif v-annotate.pl --split --cpu 8 --glsearch -s -r --nomisc --mkey sarscov2 --lowsim5seq 6 --lowsim3seq 6 --alt_fail lowscore,insertnn,deletinn --mdir /opt/vadr/vadr-models/ /data/"+items[-1]+".trimmed.fasta -f /data/"+"vadr_results", shell=True, stdout=out_log, stderr=err_log, check=True)
        subprocess.run("singularity exec -B "+"${mypath}"+"/assembly:/data docker://staphb/vadr:1.3 v-annotate.pl --split --cpu 8 --glsearch -s -r --nomisc --mkey sarscov2 --lowsim5seq 6 --lowsim3seq 6 --alt_fail lowscore,insertnn,deletinn --mdir /opt/vadr/vadr-models/ /data/"+items[-1]+".trimmed.fasta -f /data/"+"vadr_results", shell=True, stdout=out_log, stderr=err_log, check=True)

        #Parse through VADR outputs to get PASS or REVIEW flag
        vadr_flag = ''
        with open("${mypath}"+"/assembly/vadr_results/vadr_results.vadr.pass.list", 'r') as p_list:
            result = p_list.readline()
            result = result.rstrip()
            if result == items[-1]:
                vadr_flag = 'PASS'

        with open("${mypath}"+"/assembly/vadr_results/vadr_results.vadr.fail.list", 'r') as f_list:
            result = f_list.readline()
            result = result.rstrip()
            if result == items[-1]:
                vadr_flag = 'REVIEW'

        #Copy VADR error report to main analysis folder for easier review
        if vadr_flag == 'REVIEW':
            subprocess.run("cp " + "${mypath}"+"/assembly/vadr_results/vadr_results.vadr.alt.list "+"${params.output}"+"/vadr_error_reports/"+items[-1]+".vadr.alt.list", shell=True, check=True)
            #subprocess.run('mv vadr_error_reports/vadr_results.vadr.alt.list vadr_error_reports/' + s + '.vadr.alt.list', shell=True, check=True)
            #print("cp " + "${mypath}"+"/assembly/vadr_results/vadr_results.vadr.alt.list "+"${params.output}"+"/vadr_error_reports/"+items[-1]+".vadr.alt.list")

        #Run pangolin
        ###--pango_path /apps/staphb-toolkit/containers/pangolin_4.1.2-pdata-1.13.sif --pangolin v4.1.2 --pangolin_data v1.13
        pangolin = "v4.1.2_pdata-v1.13"
        #subprocess.run('singularity exec -B '+"${mypath}"+"/assembly:/data /apps/staphb-toolkit/containers/pangolin_4.1.2-pdata-1.13.sif" + ' pangolin -o /data /data/'+items[-1]+".consensus.fa", shell=True, check=True)
        subprocess.run('singularity exec -B '+"${mypath}"+"/assembly:/data docker://staphb/pangolin:4.1.2-pdata-1.13" + ' pangolin -o /data /data/'+items[-1]+".consensus.fa", shell=True, check=True)

        #Get lineage
        proc = subprocess.run("tail -n 1 "+"${mypath}"+"/assembly/lineage_report.csv | cut -d \',\' -f 2", shell=True, check=True, capture_output=True, text=True)
        lineage = proc.stdout.rstrip()
        #print(lineage)

        #Run nextclade
        #subprocess.run("singularity exec -B "+"${mypath}"+"/assembly:/data /apps/staphb-toolkit/containers/nextclade_2021-03-15.sif nextclade --input-fasta /data/"+items[-1]+".consensus.fa --output-csv /data/nextclade_report.csv", shell=True, check=True)
        subprocess.run("singularity exec -B "+"${mypath}"+"/assembly:/data docker://nextstrain/nextclade:latest nextclade dataset get --name sars-cov-2 --output-dir /data/dataset", shell=True, check=True)
        subprocess.run("singularity exec -B "+"${mypath}"+"/assembly:/data docker://nextstrain/nextclade:latest nextclade run --input-dataset /data/dataset --output-csv /data/nextclade_report.csv /data/"+items[-1]+".consensus.fa", shell=True, check=True)

        #Parse nextclade output and screen for sotc
        sotc_v = "${params.sotc}".split(',')
        with open("${mypath}"+"/assembly/nextclade_report.csv", 'r') as nc:
            header = nc.readline()
            c_results = nc.readline()
            c_results = c_results.rstrip()
            data = c_results.split(',')
            #print(data)
            sotc = []
            for v in sotc_v:
                print(v)
                if v in data:
                    sotc.append(v)
            sotc_out = (',').join(sotc)
        out_log.close()
        err_log.close()
    else:
        vadr_flag = 'NA'
        lineage = 'NA'
        sotc_out = 'NA'
        
    with open("${mypath}"+"/report.txt", 'w') as report:
        header = ['sampleID', 'reference', 'start', 'end', 'num_raw_reads', 'num_clean_reads', 'num_mapped_reads', 'percent_mapped_clean_reads', 'cov_bases_mapped', 'percent_genome_cov_map', 'mean_depth', 'mean_base_qual', 'mean_map_qual', 'assembly_length', 'numN', 'percent_ref_genome_cov', 'VADR_flag', 'QC_flag', 'pangolin_version', 'lineage', 'SOTC']
        report.write('\t'.join(map(str,header)) + '\n')
        results = [items[-1], ref_name, start, end, raw_reads, clean_reads, reads_mapped, percent_map, cov_bases, cov, depth, baseq, mapq, num_bases, ns, pg, vadr_flag, qc_flag, pangolin, lineage, sotc_out]
        report.write('\t'.join(map(str,results)) + '\n')
    /$
}
