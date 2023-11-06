for filepath in Result_*/CRISPRFinderProperties/*/Spacers/*; do

	# echo "$filepath"
	sample="${filepath#Result_*}" # Removes the "Result" part from the name
	sample="${sample%%/CRISPR*}"  # Removes everything including/after CRISPR, just keeps the sample name

	# echo "$sample"
	sequence="${filepath#${filepath%k99*}}" # Removes everything until k99
	sequence="${sequence%%_properties*}" # Removes everything from properties onward

	# echo "$sequence"
	alignment="${filepath##Result*/}" # this just keeps the DR/Spacer portion of filepath

	# echo "$alignment"
	content="$(cat $filepath)"

	# echo ""$content""

	# Send the information in the appropriate order to the csv file
	echo "$sample,$sequence,$alignment,"$content",$filepath" >> results_DRs.csv
done
