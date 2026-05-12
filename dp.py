#!/usr/bin/env python3

from os import getenv
from time import sleep
from dotenv import load_dotenv

import requests

load_dotenv()

# load the API key from the .env file, or use a placeholder value if the .env file is not found or the API_KEY variable is not set
apiKey = getenv("cyalq0qv5pmrsg9eckbkhq0lduofuu") # replace this with your actual phamerator API key
if apiKey == "cyalq0qv5pmrsg9eckbkhq0lduofuu":
apiKey = getenv("API_KEY") or "your phamerator API key here"
if apiKey == "your phamerator API key here":
    print("Error: API key not found. Create an account on https://phamerator.org and generate an API key in the account settings, then add it to a .env file in the same directory as this script with the following line:\n" \
    "API_KEY=your phamerator API key here")
    print("Please replace the apiKey variable with your actual phamerator API key,\n" \
    "which can be obtained by creating an account on https://phamerator.org and generating\n" \
    "a key in the account settings.")
    exit()

dataset = "Actino_Draft"

subcluster = "K2" # use the subcluster that you want to analyze, for example "G2" for the G2 subcluster of actinobacteriophages

match = 3
mismatch = -2
gap = -5

headers = {
    "Authorization": f"Bearer {apiKey}",
    "Accept": "application/json"
}
base_url = "https://phamerator.org/api/"

def getPhameratorData(endpoint, args=""):
  sleep(1) # add a delay of 1 second between API calls to avoid hitting the rate limit
  url = base_url + dataset + "/" + endpoint + args
  return requests.get(url, headers=headers).json()

# for a pair of phages, create a gene matrix that contains the phamily names of the genes for each phage, and then create a score matrix and an arrow matrix to store the scores and directions for the dynamic programming algorithm
def createDPMatrix(phage1, phage2):
    geneMatrix = []
    scoreMatrix = []
    arrowMatrix = []

    for phage in [phage1, phage2]:
        genes = getPhameratorData("genes/", phage["phagename"])
        # for each gene in the list of genes for this phage, append it's phamily name to the geneMatrix
        geneMatrix.append([gene["phamName"] for gene in genes])

    # create an empty score matrix with a row for each gene in geneMatrix[0] and a column for each gene in geneMatrix[1] plus one extra column and one extra row for the base case of the dynamic programming algorithm
    # the score matrix should be initialized with None values, except for the first row and first column which should be initialized with 0 values
    scoreMatrix = [[None for _ in range(len(geneMatrix[1]) + 1)] for _ in range(len(geneMatrix[0]) + 1)]

    # Initialize the first row and first column with 0 values
    for i in range(len(scoreMatrix)):
        scoreMatrix[i][0] = i * gap
    for j in range(len(scoreMatrix[0])):
        scoreMatrix[0][j] = j * gap

    # for each cell in the score matrix, calculate the score based on the following rules:
    # if the gene in the first phage is the same as the gene in the second phage, the score is 1 plus the score of the cell diagonally above and to the left
    # if the gene in the first phage is different from the gene in the second phage, the score is the maximum of the cell diagonally above and to the left, the score of the cell directly above -1 gap penalty, and the score of the cell directly to the left -1 gap penalty
    # the arrow matrix should be updated to indicate the direction of the score (diagonal, up, or left) for each cell in the score matrix
    arrowMatrix = [[None for _ in range(len(geneMatrix[1]) + 1)] for _ in range(len(geneMatrix[0]) + 1)]
    # for each row...
    for i in range(1, len(scoreMatrix)):
        # for each column...
        for j in range(1, len(scoreMatrix[0])):
            # if the gene phamilies in the two phages match, the score is match + the score of the cell diagonally above and to the left, and the arrow is "D" for diagonal
            if geneMatrix[0][i-1] == geneMatrix[1][j-1]:
                scoreMatrix[i][j] = scoreMatrix[i-1][j-1] + match # match 
                arrowMatrix[i][j] = "D"
            else:
                scores = [
                    scoreMatrix[i-1][j-1] + mismatch,  # diagonal: mismatch
                    scoreMatrix[i-1][j] + gap,     # up: gap in sequence 2
                    scoreMatrix[i][j-1] + gap      # left: gap in sequence 1
                ]
                max_score = max(scores)
                scoreMatrix[i][j] = max_score
                arrowMatrix[i][j] = ["D", "U", "L"][scores.index(max_score)]
    
    for row in scoreMatrix:
        print("\t".join(f"{score:>5}" for score in row))

    # backtrack through the score matrix and arrow matrix to align the two lists of phamilies, starting from the bottom right cell and moving to the top left cell
    alignedPhams1 = []
    alignedPhams2 = []
    i = len(scoreMatrix) - 1
    j = len(scoreMatrix[0]) - 1
    while i > 0 or j > 0:
        if arrowMatrix[i][j] == "D":
            alignedPhams1.append(geneMatrix[0][i-1])
            alignedPhams2.append(geneMatrix[1][j-1])
            i -= 1
            j -= 1
        elif arrowMatrix[i][j] == "U":
            alignedPhams1.append(geneMatrix[0][i-1])
            alignedPhams2.append("---")
            i -= 1
        else:
            alignedPhams1.append("---")
            alignedPhams2.append(geneMatrix[1][j-1])
            j -= 1

    # print the aligned lists of phamilies in the correct order (reverse the lists before printing)
    # print two tab-separated lines, one for each aligned list of phamilies, with the phamily names in the same column for both lines
    # align the columns using tabs to separate the phamily names and gap characters and setting the width of each column to be the same for both lines
    alignedPhams1.reverse()
    alignedPhams2.reverse()
    max_width = max(len(str(pham)) for pham in alignedPhams1 + alignedPhams2)
    print("\t".join(f"{pham:<{max_width}}" for pham in alignedPhams1))
    print("\t".join(f"{pham:<{max_width}}" for pham in alignedPhams2))

    return {phage1["phagename"]: alignedPhams1, phage2["phagename"]: alignedPhams2}

subclusters = ["P5", "AC"]

for subcluster in subclusters:
    # for each possible pair of phages in the specified subcluster, get the list of genes for each phage and use a dynamic programming algorithm to align the two lists of genes based on their phamily names, with a score of 1 for a match and a gap penalty of -1 for a mismatch or gap. The output should be two aligned lists of phamilies, with gaps represented by "---", and the columns should be aligned using tabs to separate the phamily names and gap characters and setting the width of each column to be the same for both lines.
    phages = [phage for phage in getPhameratorData("phagesbysubcluster/", subcluster)]

    alignments=[]
    # put all the phages in the specified subcluster into a list and then iterate through all possible pairs of phages in the list to compare them using the dynamic programming algorithm
    for i in range(len(phages)):
        for j in range(i + 1, len(phages)):
            # print(f"Comparing {phages[i]['phagename']} and {phages[j]['phagename']}...")
            dpMatrix = createDPMatrix(phages[i], phages[j])
            alignments.append((dpMatrix[phages[i]["phagename"]], dpMatrix[phages[j]["phagename"]]))

    # save alignments to a file named "alignments_{subcluster}.tsv" in the same directory as this script, with the aligned lists of phamilies for each pair of phages printed in the same format as the output of the createDPMatrix function, use the match, mismatch and gap scores to show the score of each aligned column, and also show the total score of each alignment at the end of the aligned lists of phamilies, and separate each pair of aligned lists of phamilies with a blank line for readability
    with open(f"alignments_{subcluster}.tsv", "w") as f:
        for alignment in alignments:
            phams1, phams2 = alignment
            # begin each alignment line with the names of the two phages that were aligned, separated by a tab, and then print the aligned lists of phamilies in the same format as the output of the createDPMatrix function, using tabs to separate the phamily names and gap characters and setting the width of each column to be the same for both lines
            max_width = max(len(str(pham)) for pham in phams1 + phams2)
            f.write(f"{phages[alignments.index(alignment) // (len(phages) - 1)]['phagename']}\t{phages[alignments.index(alignment) % (len(phages) - 1) + 1]['phagename']}\n")
            f.write("\t".join(f"{pham:<{max_width}}" for pham in phams1) + "\n")
            f.write("\t".join(f"{pham:<{max_width}}" for pham in phams2) + "\n")
            # write the total score of the alignment, using match, mismatch and gap scores to calculate the score of each aligned column and the total score of the alignment
            column_scores = []
            for p1, p2 in zip(phams1, phams2):
                if p1 == p2 and p1 != "---":
                    column_scores.append(match)  # match
                elif p1 == "---" or p2 == "---":
                    column_scores.append(gap)  # gap penalty
                else:
                    column_scores.append(mismatch)  # mismatch penalty
            f.write("\t".join(f"{score:<{max_width}}" for score in column_scores) + "\n")
            
            # total_score is the sum of the column scores
            total_score = sum(column_scores)

            f.write(f"Total Score: {total_score}\n")
            f.write("\n")  # separate each pair of aligned lists with a blank line
        # close the file after writing all the alignments
    f.close()
    ungappedAlignments = []

    for alignment in alignments:
        phams1, phams2 = alignment
        for index in range(len(phams1)):
            if phams1[index] != "---" and phams2[index] != "---":
                ungappedAlignments.append((phams1[index], phams2[index]))


    totalCounts = {}
    for q, s in ungappedAlignments:
        # count the number of times q or s appears in the ungapped alignments
        if q not in totalCounts:
            totalCounts[q] = 0
        if s not in totalCounts:
            totalCounts[s] = 0
        totalCounts[q] += 1
        totalCounts[s] += 1

    print (totalCounts)

    # for each phamily, count the number of times it is aligned with each of the phamilies in the ungapped alignments and store the counts in a dictionary of dictionaries, where the keys of the outer dictionary are the phamily names and the values are dictionaries with the keys being the phamily names that they are aligned with and the values being the counts of how many times they are aligned together
    alignmentCounts = {}
    for q, s in ungappedAlignments:
        if q not in alignmentCounts:
            alignmentCounts[q] = {}
        if s not in alignmentCounts[q]:
            alignmentCounts[q][s] = 0
        alignmentCounts[q][s] += 1

        if s not in alignmentCounts:
            alignmentCounts[s] = {}
        if q not in alignmentCounts[s]:
            alignmentCounts[s][q] = 0
        alignmentCounts[s][q] += 1

    print (alignmentCounts)
