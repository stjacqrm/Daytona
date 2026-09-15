process pystats2 {
    input:
        //val mypath
        val x
    output:
        stdout
        //val mypath
        //path "pyoutputs.txt", emit: pyoutputs        
    $/
    #!/usr/bin/env python3
    import subprocess
    from Bio import SeqIO
    
    subprocess.run("mkdir -p " + "${params.output}"+"/"+"${x}", shell=True, check=True)
    #subprocess.run('mkdir -p vadr_error_reports_clearlabs', shell=True, check=True)
    #subprocess.run('mkdir -p '+"${params.output}"+"/"+"${x}"+"/"+'assemblies_pass', shell=True, check=True) 
    
    out_log = open("${params.output}"+"/"+"${x}"+"/" + "${x}" + ".out", 'w')
    err_log = open("${params.output}"+"/"+"${x}"+"/" + "${x}"+ ".err", 'w')
    
    #Get number of raw reads
    proc_r = subprocess.run('cat ' + "${params.input}"+"/"+"${x}" + "*.fastq | wc -l", shell=True, capture_output=True, text=True, check=True)
    wc_out_1 = proc_r.stdout.rstrip()
    clean_reads = int(wc_out_1) / 4
    #print(clean_reads)
    
    #Get bam file name
    proc_b = subprocess.run("ls " + "${params.bams}"+"/"+"${x}" + "*.bam", shell=True, capture_output=True, text=True, check=True)
    bam_file = proc_b.stdout.rstrip()
    items = bam_file.strip().split("/")
    #print(items[-1])
    
    #Generate bam with only mapped reads
    subprocess.run('singularity exec -B '+"${params.bams}"+':/data ' + "${params.samtools_docker}" + ' samtools view -F 4 -u -h /data/' + items[-1] + ' | singularity exec -B '+"${params.bams}"+':/data ' + "${params.samtools_docker}" + ' samtools sort > ' + "${params.output}"+"/"+"${x}"+"/"+"${x}"+ '.mapped.sorted.bam', shell=True, check=True)
    subprocess.run('singularity exec -B '+"${params.output}"+"/"+"${x}"+':/data ' + "${params.samtools_docker}" + ' samtools index /data/' + "${x}"+'.mapped.sorted.bam', shell=True, check=True)    
    
    #Run samtools coverage to get map stats
    subprocess.run('singularity exec -B '+"${params.output}"+"/"+"${x}"+':/data ' + "${params.samtools_docker}" + ' samtools coverage /data/' + "${x}"+ '.mapped.sorted.bam -o /data/' + "${x}" + '.coverage.txt', shell=True, stdout=out_log, stderr=err_log, check=True)
    
    #Get map stats
    with open("${params.output}"+"/"+"${x}"+"/"+"${x}"+'.coverage.txt', 'r') as cov_report:
        header = cov_report.readline()
        header = header.rstrip()
        stats = cov_report.readline()
        stats = stats.rstrip()
        stats = stats.split()
        ref_name = stats[0]
        start = stats[1]
        end = stats[2]
        reads_mapped = stats[3]
        cov_bases = stats[4]
        cov = stats[5]
        depth = stats[6]
        baseq = stats[7]
        mapq = stats[8]
    
    #Get percentage of mapped reads/reads
    percent_map = "%0.4f"%(((int(reads_mapped)/int(clean_reads)))*100)
    
    #Get fasta file name
    proc_f = subprocess.run('ls ' + "${params.assemblies}"+"/"+ "${x}" + '*.fasta', shell=True, capture_output=True, text=True, check=True)
    fasta_file = proc_f.stdout.rstrip()
    
    #Gather QC metrics for consensus assembly
    for record in SeqIO.parse(fasta_file, "fasta"):
        seq = record.seq
        num_bases = len(seq)
        ns = seq.lower().count('n')
        called = num_bases - ns
        pg = "%0.4f"%((called/int(end))*100)
    #print(pg)
    
    #Copy assembly to sample dir in order to rename and run VADR
    subprocess.run('cp ' + "${params.assemblies}"+"/"+ "${x}"+ '*.fasta ' + "${params.output}"+"/"+"${x}"+"/" + "${x}" + '.consensus.fa', shell=True, check=True)
    #Rename header in fasta to just sample name
    subprocess.run('sed -i \'s/^>.*/>' + "${x}" + '/\' ' + "${params.output}"+"/"+"${x}"+"/" + "${x}" + '.consensus.fa', shell=True, check=True)
    
    #QC flag
    pg_flag = ''
    dp_flag = ''
    qc_flag = ''
    if float(pg) < 89.5:
        pg_flag = 'FAIL: Percent genome < 90%'
        qc_flag = qc_flag + pg_flag
    else:
        if float(depth) < 100:
            dp_flag = 'FAIL: Mean read depth < 100x'
            qc_flag = qc_flag + dp_flag
        if qc_flag == '':
            qc_flag = qc_flag + 'PASS'
    print(qc_flag)

    #Copy all assemblies
    subprocess.run('mkdir -p '+"${params.output}"+"/"+'assemblies_total', shell=True, check=True)
    subprocess.run('cp ' + "${params.output}"+"/"+"${x}"+"/" + "${x}" + '.consensus.fa '+"${params.output}"+"/"+'assemblies_total', shell=True, check=True)

    #Run pangolin
    #subprocess.run('singularity exec -B '+"${params.output}"+"/"+"${x}"+':/data ' + "${params.pangolin_docker}" + ' pangolin /data/' + "${x}" + '.consensus.fa -o '+"${params.output}"+"/"+"${x}", shell=True, check=True)
    subprocess.run('singularity exec -B '+"${params.output}"+"/"+"${x}"+':/data ' + "${params.pangolin_docker}" + ' pangolin /data/' + "${x}" + '.consensus.fa -o /data', shell=True, check=True)

    #Get lineage
    proc = subprocess.run('tail -n 1 '+"${params.output}"+"/"+"${x}"+'/lineage_report.csv | cut -d \',\' -f 2', shell=True, check=True, capture_output=True, text=True)
    lineage = proc.stdout.rstrip()
    print(lineage)
    
    #Copy passing assemblies to assemblies folder and vcf files to variants folder
    if qc_flag == 'PASS':
        subprocess.run('cp ' + "${params.output}"+"/"+"${x}"+"/" + "${x}" + '.consensus.fa '+"${params.output}"+"/"+'assemblies_pass/', shell=True, check=True)
        
        #Run VADR
        subprocess.run('singularity exec -B '+"${params.output}"+"/"+"${x}"+':/data ' + "${params.vadr_docker}" + ' /opt/vadr/vadr/miniscripts/fasta-trim-terminal-ambigs.pl --minlen 50 --maxlen 30000 /data/' + "${x}" + '.consensus.fa > '+"${params.output}"+"/"+"${x}"+"/"+"${x}"+'.trimmed.fasta', shell=True, stdout=out_log, stderr=err_log, check=True)
        subprocess.run('singularity exec -B '+"${params.output}"+"/"+"${x}"+':/data ' + "${params.vadr_docker}" + ' v-annotate.pl --noseqnamemax --split --cpu 8 --glsearch -s -r --nomisc --mkey sarscov2 --lowsim5seq 6 --lowsim3seq 6 --alt_fail lowscore,insertnn,deletinn --mdir /opt/vadr/vadr-models/ /data/' + "${x}" + '.trimmed.fasta -f /data/vadr_results', shell=True, stdout=out_log, stderr=err_log, check=True)

        #Parse through VADR outputs to get PASS or REVIEW flag
        vadr_flag = ''
        with open("${params.output}"+"/"+"${x}"+"/"+ 'vadr_results/vadr_results.vadr.pass.list', 'r') as p_list:
            result = p_list.readline()
            result = result.rstrip()
            if result == "${x}":
                vadr_flag = 'PASS'
                
        with open("${params.output}"+"/"+"${x}"+"/" + 'vadr_results/vadr_results.vadr.fail.list', 'r') as f_list:
            result = f_list.readline()
            result = result.rstrip()
            if result == "${x}":
                vadr_flag = 'REVIEW'
        #print(vadr_flag)
        
        #Copy VADR error report to main analysis folder for easier review
        if vadr_flag == 'REVIEW':
            subprocess.run('cp ' + "${params.output}"+"/"+"${x}"+"/" + 'vadr_results/vadr_results.vadr.alt.list '+ "${params.output}"+"/"+'vadr_error_reports_clearlabs/'+ "${x}" + '.vadr.alt.list', shell=True, check=True)
            #subprocess.run('mv '+ "${params.output}"+"/"+'vadr_error_reports_clearlabs/vadr_results.vadr.alt.list '+"${params.output}"+"/"+'vadr_error_reports_clearlabs/' + "${x}" + '.vadr.alt.list', shell=True, check=True)


        #Run nextclade
        subprocess.run('singularity exec -B '+"${params.output}"+"/"+"${x}"+':/data ' + "${params.nextclade_docker}" + ' nextclade --input-fasta /data/' + "${x}" + '.consensus.fa --output-csv /data/' + 'nextclade_report.csv', shell=True, check=True)
        
        #Parse nextclade output and screen for sotc
        sotc_v = "${params.sotc}".split(',')
        with open("${params.output}"+"/"+"${x}" + '/nextclade_report.csv', 'r') as nc:
            header = nc.readline()
            c_results = nc.readline()
            c_results = c_results.rstrip()
            data = c_results.split(',')
            sotc = []
            for v in sotc_v:
                if v in data:
                    sotc.append(v)
            sotc_out = (',').join(sotc)
    else:
        vadr_flag = 'NA'
        #lineage = 'NA'
        sotc_out = 'NA'
    
    pangolin = "${params.pangolin_version}"+'_pdata-'+"${params.pangolin_data_version}"
    with open("${params.output}"+"/"+"${x}"+"/report.txt", 'w') as report:
        header = ['sampleID', 'reference', 'start', 'end', 'num_clean_reads', 'num_mapped_reads', 'percent_mapped_clean_reads', 'cov_bases_mapped', 'percent_genome_cov_map', 'mean_depth', 'mean_base_qual', 'mean_map_qual', 'assembly_length', 'numN', 'percent_ref_genome_cov', 'VADR_flag', 'QC_flag', 'pangolin_version', 'lineage', 'SOTC']
        report.write('\t'.join(map(str,header)) + '\n')
        results = ["${x}", ref_name, start, end, int(clean_reads), reads_mapped, percent_map, cov_bases, cov, depth, baseq, mapq, num_bases, ns, pg, vadr_flag, qc_flag, pangolin, lineage, sotc_out]
        report.write('\t'.join(map(str,results)) + '\n')
    /$
}
