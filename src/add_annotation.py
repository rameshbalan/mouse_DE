#!/usr/bin/python3

#Declaring list and dictionary
transcript_list = []
transcript_ann = {}

#Reading the file with has significant differential expression:
with open("significant_transcripts.txt","r") as transcript_file:

	for line in transcript_file:

		transcript_list.append(line)

#Adding annotation based on the gff file
with open("ref_GRCm38.p4_top_level.gff3","r") as gff_file:

	for line in gff_file:

		if not line.startswith("#"):

			split_line = line.split("\t")

			if split_line[2] == "mRNA":

				split_col9 = split_line[8].split(";")

				transcript_id = split_col9[-1].split("=")

				if transcript_id[1] in transcript_list:

					transcript_ann[transcript_id[1].rstrip()] = split_col9[-2].split("=")[1]

print(transcript_ann)

#Writing the transcript ID and annotation to a new file
with open("significant_transcripts_annotated.txt","w") as outfile:

	for key, value in transcript_ann.items():

		outfile.write("%s\t%s\n" % (key,value))