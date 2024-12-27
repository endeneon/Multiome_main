for eachfile in *.bam
do
	echo $eachfile
	samtools index -@ 20 $eachfile
done

