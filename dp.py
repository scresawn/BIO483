#!/usr/bin/env python3

from os import getenv
from dotenv import load_dotenv

import requests

load_dotenv()

# load the API key from the .env file, or use a placeholder value if the .env file is not found or the API_KEY variable is not set
apiKey = getenv("cyalq0qv5pmrsg9eckbkhq0lduofuu") # replace this with your actual phamerator API key
if apiKey == "cyalq0qv5pmrsg9eckbkhq0lduofuu":
    print("Please replace the apiKey variable with your actual phamerator API key,\n" \
    "which can be obtained by creating an account on https://phamerator.org and generating\n" \
    "a key in the account settings.")
    exit()

dataset = "Actino_Draft"

subcluster = "K2" # use the subcluster that you want to analyze, for example "G2" for the G2 subcluster of actinobacteriophages

headers = {
    "Authorization": f"Bearer {apiKey}",
    "Accept": "application/json"
}
base_url = "https://phamerator.org/api/"

def getPhameratorData(endpoint, args=""):
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
        scoreMatrix[i][0] = 0
    for j in range(len(scoreMatrix[0])):
        scoreMatrix[0][j] = 0

    # for each cell in the score matrix, calculate the score based on the following rules:
    # if the gene in the first phage is the same as the gene in the second phage, the score is 1 plus the score of the cell diagonally above and to the left
    # if the gene in the first phage is different from the gene in the second phage, the score is the maximum of the cell diagonally above and to the left, the score of the cell directly above -1 gap penalty, and the score of the cell directly to the left -1 gap penalty
    # the arrow matrix should be updated to indicate the direction of the score (diagonal, up, or left) for each cell in the score matrix
    arrowMatrix = [[None for _ in range(len(geneMatrix[1]) + 1)] for _ in range(len(geneMatrix[0]) + 1)]
    for i in range(1, len(scoreMatrix)):
        for j in range(1, len(scoreMatrix[0])):
            if geneMatrix[0][i-1] == geneMatrix[1][j-1]:
                scoreMatrix[i][j] = scoreMatrix[i-1][j-1] + 1
                arrowMatrix[i][j] = "D"
            else:
                scores = [
                    scoreMatrix[i-1][j-1],  # diagonal
                    scoreMatrix[i-1][j] - 1,     # up
                    scoreMatrix[i][j-1] - 1      # left
                ]
                max_score = max(scores)
                scoreMatrix[i][j] = max_score
                arrowMatrix[i][j] = ["D", "U", "L"][scores.index(max_score)]

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

# for each possible pair of phages in the specified subcluster, get the list of genes for each phage and use a dynamic programming algorithm to align the two lists of genes based on their phamily names, with a score of 1 for a match and a gap penalty of -1 for a mismatch or gap. The output should be two aligned lists of phamilies, with gaps represented by "---", and the columns should be aligned using tabs to separate the phamily names and gap characters and setting the width of each column to be the same for both lines.
phages = [phage for phage in getPhameratorData("phagesbysubcluster/", subcluster)]

# put all the phages in the specified subcluster into a list and then iterate through all possible pairs of phages in the list to compare them using the dynamic programming algorithm
for i in range(len(phages)):
    for j in range(i + 1, len(phages)):
        print(f"Comparing {phages[i]['phagename']} and {phages[j]['phagename']}...")
        createDPMatrix(phages[i], phages[j])