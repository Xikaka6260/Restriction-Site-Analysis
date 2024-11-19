import sys

file1=open(sys.argv[1])
file2=open(sys.argv[2])
#argv 1 should be the fasta file
#argv 2 should be the enzyme file


file1_lines=file1.readlines()
file2_lines=file2.readlines()
#this step turns every line read from both files into str format

file1.close()
file2.close()

sequences_list=[line.rstrip() for line in file1_lines]	#removing newline characters
general_enzymes_list=[line.rstrip().replace("^","") for line in file2_lines]
#print(sequences_list)
#print(general_enzymes_list)

#Making a new variable that contain only the sequences without the header
actual_sequence=sequences_list[1:]
#print(sequences)

print()
print("Restriction enzyme analysis of sequence from file " + sys.argv[1] + ".\n" + "Cutting with enzymes found in file " + sys.argv[2] + ".")
print("---"*43)
header=sequences_list[0]
print("Sequence name :" + header.replace(">",""))
presequence=",".join(actual_sequence)	#the sequence were originally separated,here joining back together
sequence=presequence.replace(",","")	#removing all the commas
print("Sequence is " + str(len(sequence)) + " bases long.")
print("---"*43)
print()


#Making a dictionary and assignment enzymes as key, and their restriction site as value
enzymes_dict={}

for j in general_enzymes_list:
	enzyme, restriction_site = j.split(";")	#pattern of the enzyme and cutting sites are separated in ;
	enzymes_dict[enzyme]=restriction_site

#print(enzymes_dict)


#looping through all of the enzymes possibles and retrieve their corresponding restriction site
restriction_site_list=[]
for enzyme in enzymes_dict:
	restriction_site=enzymes_dict.get(enzyme)
	restriction_site_list.append(restriction_site)
#print(restriction_site_list)


#Trying when each restriction_site is extracted from dictionary previously
#This step has to be done because dictionary has no specific location for each key
def enzyme_cutting_site (sequence,restriction_site_list):
	cut_sites=[]	#creating a list to store each cut_locations for downstream use

	for w in restriction_site_list:		#loop through every restriction_sites previously stated
		start=0
		cut_locations=[]		#keep storing all cutting sites each run

		while start < len(sequence):
			position=sequence.find(w,start)
						#find the location where restriction_site matches and keep track of location
			if position != -1:	#run this loop only when restriction is found
				start = position +1
				cut_locations.append(str(position))

			else :
				break
		cut_sites.append([w,cut_locations])

	return cut_sites

		#print(f" When restriction_site is {w} cut locations are {cut_locations}")


#Now I stored all the cutting sites for each enzyme into a list, this is a list in a list situation, for downstream use
cut_sites=(enzyme_cutting_site(sequence,restriction_site_list))
#print(cut_sites)



fragment_count=len(cut_sites[0][1])+1	#base on the cut_sites list, the cutting sites for this enzyme is the second list


#Looping to find which enzyme's restriction site matches, and extract the enzyme name
enzyme1=None
for enzyme, restriction_site in enzymes_dict.items():
	if restriction_site == cut_sites[0][0]:
		enzyme1=enzyme
		break
#print(enzyme1)

#print(cut_locations)
#print(type(cut_locations[0]))



#Extract second list from the cut_sites list, these are all cutting sites only for each enzymes
extracted_cut_sites=[]
for q in cut_sites:
	extracted_cut_sites.append((q[1]))

#print(extracted_cut_sites)



#Turning all these cutting sites to integers
for r  in extracted_cut_sites:
	for u in range(len(r)):
		r[u]=int(r[u])
#print(extracted_cut_sites)


cut_locations=extracted_cut_sites[0]


#Making a function that can match through all cutting sites for enyzme1 and provide the formatted fragments
def process_fragments(sequence, cut_locations):

	starting_location = 0	# keeping count of the current location of each fragments
	fragments = []		# making an empty list to store my fragments


	#loop through every cutting site
	for cut in cut_locations:
		fragment = sequence[starting_location:cut]		#extract each fragment from the sequence
		fragments.append((fragment, starting_location, cut))	#store each fragment into my empty list
		starting_location = cut  # Update starting location


   	 # Add the last fragment when there is no more cutting sites available
	fragment = sequence[starting_location:]
	fragments.append((fragment, starting_location, len(sequence)))


	 # Print fragments and their lengths

	for b, (fragment, start, end) in enumerate(fragments):	#for every fragment in the list
		result=[]
		for t in range(0, len(fragment),60):		#for every line of the fragment
			display_fragment=fragment[t:t+60]	#show every line only contain 60 bases
			arranged_fragment=" ".join(display_fragment[i:i+10] for i in range(0, len(display_fragment),10))
								#for every line, show 10 subgroups each
			result.append(arranged_fragment)

		#print(f"Length - {len(fragment)}")
		print(f"Starting Position {start}"+ "  " + f"Fragment length - {len(fragment)}")

		for res in result:		#print out every re-formatted fragments
			print(res)
		print("\n")




#Print fragments by calling this function giving it the sequence we uses, and the cut_locations for each enzyme
if cut_sites:
	print("There are " + str(len(cut_sites[0][1])) + " cutting sites for " + str(enzyme1) + ".  Cutting at " + str(cut_sites[0][0])+ "." + "\nThere are " + str(fragment_count) + " fragments: ")
	print("\n")
	process_fragments(sequence, cut_locations)
else:
	print("No Cutting Site for " + enzyme1)




#------------------------------------------------------------Second Enzyme-------------------------------------------------------
print("---"*43)
print()
enzyme2 = None
# Loop through the enzyme diction to find the match for second enzyme
for enzyme, restriction_site in enzymes_dict.items():
	if restriction_site == cut_sites[1][0]:  # Match the restriction site
		enzyme2 = enzyme
		break  # Exit the loop after finding the match

#print(enzyme2)

cut_locations=extracted_cut_sites[1]	#Getting cut_location for second enzyme


if cut_locations:
        process_fragments(sequence, cut_locations)
else:
     	print("No Cutting Site for " + enzyme2)




#-----------------------------------------------------------Third Enzyme------------------------------------------------------
print("---"*43)
print()
enzyme3 = None
# Loop through the enzyme diction to find the match for third enzyme
for enzyme, restriction_site in enzymes_dict.items():
	if restriction_site == cut_sites[2][0]:  # Match the restriction site
		enzyme3 = enzyme
		break  # Exit the loop after finding the match

#print(enzyme3)

cut_locations=extracted_cut_sites[2]	#cutting site for this enzyme is the 3rd list

fragment_count=len(cut_locations)+1

if cut_locations:
	print("There are " + str(len(cut_locations)) + " cutting sites for " + enzyme3 + ".  Cutting at " + str(cut_sites[2][0]) + "." + "\nThere are " + str(fragment_count) + " fragments: ")
	print("\n")
	process_fragments(sequence, cut_locations)
else:
	print("No Cutting Site for " + enzyme3)



#-------------------------------------------------------------Fourth Enzyme---------------------------------------------------

print("---"*43)
print()
enzyme4 = None
# Loop through the enzyme diction to find the match for second enzyme
for enzyme, restriction_site in enzymes_dict.items():
        if restriction_site == cut_sites[3][0]:  # Match the restriction site
                enzyme4 = enzyme
                break  # Exit the loop after finding the match

#print(enzyme3)

cut_locations=extracted_cut_sites[3]    #cutting site for this enzyme is the 3rd list


fragment_count=len(cut_locations)+1


if cut_locations:
	print("There are " + str(len(cut_locations)) + " cutting sites for " + enzyme4 + ".  Cutting at " + str(cut_sites[3][0]) + "." + "\nThere are " + str(fragment_count) + " fragments: ")
	print("\n")
	process_fragments(sequence, cut_locations)
else:
	print("No Cutting Site for " + enzyme4)
