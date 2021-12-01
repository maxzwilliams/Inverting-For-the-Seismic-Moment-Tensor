"""
Program for inverting the seismic moment tensor


Almost everything in this program is written from scratch
"""


## imports 
import math ## basic math operations 
import random ## random numbers
import csv
import matplotlib.pyplot as plt ## plotting
from obspy import read
import time

def readGreensFunctions(greensPath, stationInformationPath, depth):
	"""
	Reads greens functions and station information into the program
	"""
	stationDict = dict()
	with open(stationInformationPath, 'r') as file:
		reader = csv.reader(file)
		counter = 0
		for row in reader:
			if (counter != 0):
				stationDict[row[0]] = float(row[1])
			counter = counter + 1    
		    
	components = ["ZDD", "RDD", "ZDS", "RDS", "TDS", "ZSS", "RSS", "TSS", "ZEX", "REX"]
	stationList = list(stationDict.keys())
	
	g = dict()
	for station in stationList:
		g[station] = dict()
		for c in components:
			g[station][c] = []
			
	for stationIndex in range(len(stationList)):
		for componentIndex in range(len(components)):
			fileName = greensPath + "//" + stationList[stationIndex].strip() + "."+str(depth)+".0000."+components[componentIndex]
			
			try:	
				timeSeries = read(fileName, debug_headers=True)
			except:
				raise Exception("no such file " + str(fileName) + " exists")
			g[stationList[stationIndex]][components[componentIndex]]  = list(timeSeries[0].data)
			
			
	return g, stationDict

def forwardsProblemOneStationOneTime(m, g, stationDict, stationName, timeIndex):
	"""
	solves the forwards problem for one station at one point in time
	"""
	az = stationDict[stationName] ## need to convert to radians
	az = math.pi/180 * az
	
	
	ZDD = g[stationName]["ZDD"][timeIndex]
	RDD = g[stationName]["RDD"][timeIndex]
	ZDS = g[stationName]["ZDS"][timeIndex]
	RDS = g[stationName]["RDS"][timeIndex]
	TDS = g[stationName]["TDS"][timeIndex]
	ZSS  = g[stationName]["ZSS"][timeIndex]
	RSS = g[stationName]["RSS"][timeIndex]
	TSS = g[stationName]["TSS"][timeIndex]
	ZEP = g[stationName]["ZEX"][timeIndex]
	REP = g[stationName]["REX"][timeIndex]
	

	
	uz = m[0][0]*(ZSS/2 * math.cos(2*az)- ZDD/6 + ZEP/3) + m[1][1]*(-ZSS/2*math.cos(2*az) - ZDD/6 + ZEP/3) + m[2][2]*(ZDD/3 + ZEP/3) + m[0][1]*(ZSS*math.sin(2*az)) + m[0][2]*(ZDS*math.cos(az)) + m[1][2]*(ZDS*math.sin(az))
	
	
	ur = m[0][0]*(RSS/2 * math.cos(2*az) - RDD/6 + REP/3) + m[1][1]*(-RSS/2 * math.cos(2*az) - RDD/6 + REP/3) + m[2][2] * (RDD/3 + REP/3) + m[0][1] *(RSS * math.sin(2*az)) + m[0][2] * (RDS*math.cos(az)) + m[1][2]*(RDS*math.sin(az))
	
	ut = m[0][0]*(TSS/2 * math.sin(2*az)) + m[1][1]*(-TSS/2 * math.sin(2*az)) + m[0][1] * (-TSS * math.cos(az)) + m[0][2] * (TDS * math.sin(az)) + m[1][2] * (-TDS * math.cos(az))	
	

	
	return uz, ur, ut
	
def forwardsProblemOneStation(m, g, stationDict, stationName):
	"""
	Solves the forwards problem for one station
	"""	
	uzs = []
	urs = []
	uts = []
	
	counter = 0
	while True:
		try:
			uz, ur, ut = forwardsProblemOneStationOneTime(m, g, stationDict, stationName, counter)
		except:
			break
		uzs.append(uz)
		urs.append(ur)
		uts.append(ut)
		counter += 1
	return [uzs, urs, uts]
	
	
def forwardsProblem(m, g, stationDict):
	"""
	solves the forwards problem
	"""
	stations = list(stationDict.keys())
	rtn = []
	for station in stations:
		rtn.append(forwardsProblemOneStation(m, g, stationDict, station))
		
	return rtn

	
def varience(ds, dPrimes):
	"""
	Finds the varience beteen two sets of seismograms ds and dPrimes
	"""
	s = 0
	for stationIndex in range(len(ds)):
	
		for directionIndex in range(len(ds[stationIndex])):
		
			
			partialSum = 0
			
			for index in range(len(ds[stationIndex][directionIndex])):
				di = ds[stationIndex][directionIndex][index]
				diPrime = dPrimes[stationIndex][directionIndex][index]
				try:
					partialSum += math.sqrt((di - diPrime)**2)/math.sqrt(di**2)
				except:
					partialSum += math.sqrt((di - diPrime)**2)/math.sqrt(10**(-10))
			partialSum = 1/partialSum
			s += partialSum
	return s


def randomMT(momentTensorScale):
	"""
	Generates a random seismic moment tensor that respects symmetries 
	of the siesmic moment tensor
	"""
	MT = [[0,0,0],[0,0,0],[0,0,0]]
	for i in range(3):
		for j in range(i+1):
			MT[i][j] = random.uniform(-momentTensorScale, momentTensorScale)
			MT[j][i] = MT[i][j]
			
			
	return MT 

def generatePopulation(momentTensorScale, populationSize, target):
	"""
	Generates a list of random seismic moment tensors.
	"""
	population = []
	for index in range(populationSize):
		population.append(randomMT(momentTensorScale))	
	return population
	
def scorePopulation(population, g, stationDict, realData):
	"""
	Scores population for the generati
	"""
	scoresAndPop = []

	score = 0
	for individual in population:
		syntheticData = forwardsProblem(individual, g, stationDict )
		score = varience(realData, syntheticData)
		scoresAndPop.append(  [score, individual])
	return scoresAndPop
	
def breed(mother, father, noise):
	"""
	Takes two seismic moment tensors, the mother and father and a scale for random noise. 
	A child moment tensor, which is the average of the mother and father (with some added random noise) 
	is returned.
	
	"""
	child = randomMT(0)
	for i in range(3):
		for j in range(3):
			child[i][j] = (mother[i][j]+father[i][j])/2 + (mother[i][j]-father[i][j])/2 * random.uniform(-noise, noise)
	
	
	
	for i in range(3):
		for j in range(i+1):
			child[i][j] = child[j][i]
	return child
	
def selectFromPop(selector, scoredPop):
	"""
	Takes a selector (random number between 0 and 1) and a score sorted scored population (like the output of scorePopulation)
	and returns the first moment tensor in the scoredPop with a score greater than the selector. If not such 
	moment tensor exists, something has gone wrong elsewhere in the program
	
	This function is simply used to pick a moment tensor from the population.
	"""
	
	cumulativeScoredPop = []
	cumScore = 0
	for el in scoredPop:
		cumScore += el[0]
		cumulativeScoredPop.append([cumScore, el[1]])
		
	for index in range(len(cumulativeScoredPop)):
		cumulativeScoredPop[index][0] = cumulativeScoredPop[index][0]/cumScore
		
		
		
		
	counter = 0
	for el in cumulativeScoredPop:
		if (selector <= el[0]):
			return el[1]
		counter += 1
			
	raise Exception("no selection made!")
	
def breedPopulation(scoredPop, noise):
	"""
	Takes a scored population and a random noise parameter and breeds the scored population together to generate
	a new population of the same size.
	"""
	newPop = []
	for index in range(len(scoredPop)):
		selector1 = random.uniform(0, 1)
		selector2 = random.uniform(0, 1)
		
		newChild = breed(selectFromPop(selector1, scoredPop), selectFromPop(selector2, scoredPop), noise)
		newPop.append(newChild)
		
		
	return newPop
	
	
def runGA(targetMT, momentTensorScale, popSize, epochs, geneticNoise, g, stationDict):
	"""
	Runs the genetic algorithm and plots results as it goes
	"""

	population = generatePopulation(momentTensorScale, popSize, targetMT)
	averageScores = []
	scoreDeviations = []
	
	averageMTSimilaritys = [] ## similarity between the average MT and the target
	deviationofMTs = [] ## how much deviation is there in the deviation metric
	
	realData = forwardsProblem(targetMT, g, stationDict)

	for epoch in range(epochs):
		averageScore = 0
		scoredPop = scorePopulation(population, g, stationDict, realData)
	
		maxScore = -999999999999999999
		minScore = 999999999999999999
		for element in scoredPop:
			averageScore += element[0]
			
			if (element[0] > maxScore):
				maxScore = element[0]
			if (element[0] < minScore):
				minScore = element[0]
				
		scoreDeviation = maxScore-minScore
			
		averageScore = averageScore/len(scoredPop)
		upperBound = [[-1000000, -1000000, -1000000] for _ in range(3)]
		lowerBound = [[1000000, 1000000, 1000000] for _ in range(3)]
		
		averageMT = randomMT(0)
		for individual in population:
			for i in range(3):
				for j in range(3):
					averageMT[i][j] += individual[i][j]
					upperBound[i][j] = max(upperBound[i][j], individual[i][j])
					lowerBound[i][j] = min(lowerBound[i][j], individual[i][j])

					
		for i in range(3):
			for j in range(3):
				averageMT[i][j] = averageMT[i][j]/len(population)
		deviation = [[ abs(upperBound[i][j] - lowerBound[i][j]) for j in range(3)] for i in range(3)]
		
		deviationScalar = 0
		MTsimilarity = 0
		for i in range(3):
			for j in range(3):
				deviationScalar += abs(deviation[i][j])
				MTsimilarity += abs(averageMT[i][j] - targetMT[i][j])
				
		print(MTsimilarity)
		
		
				
		
		averageScores.append(averageScore)
		scoreDeviations.append(scoreDeviation)
		averageMTSimilaritys.append(MTsimilarity) ## similarity between the average MT and the target
		deviationofMTs.append(deviationScalar) ## how much deviation is there in the deviation metric
		
		
		fig, axes = plt.subplots(1, 1)
		axes.plot(averageScores)
		axes.set_xlabel("epoch")
		axes.set_ylabel("Average Score")
		plt.savefig("scores.png")
		
		fig, axes = plt.subplots(1, 1)
		axes.plot(scoreDeviations)
		axes.set_xlabel("epoch")
		axes.set_ylabel("Deviation in score")
		plt.savefig("scoreDeviations.png")
		
		fig, axes = plt.subplots(1, 1)
		axes.plot(averageMTSimilaritys)
		axes.set_xlabel("epoch")
		axes.set_ylabel("Difference between moment tensor and target")
		plt.savefig("MTsimilarity.png")
		
		fig, axes = plt.subplots(1, 1)
		axes.plot(deviationofMTs)
		axes.set_xlabel("epoch")
		axes.set_ylabel("Deviation in moment tensors")
		plt.savefig("MTDeviation.png")
		
		
		print("==========================")
		print("epoch #:", epoch)
		print("Average Moment tensor")
		
		for index in range(3):
			print(averageMT[index])
		print("target Moment tensor")
		for index in range(3):
			print(targetMT[index])
		print("Deviation of Moment tensor in population")
		for index in range(3):
			print(deviation[index])
		print("averageScore")
		print(averageScore)
		print("==========================")
		
		predictedData = forwardsProblem(averageMT, g, stationDict)	
		plotSeismograms(realData, predictedData, epoch, stationDict)
		population = breedPopulation(scoredPop, geneticNoise)

	print("Optimized Moment Tensor")
	return averageMT

	
	
def randomList(length):
	"""
	
	Generates a random list of size length. Elements are between -1 and 1.
	"""
	rtn = []
	for index in range(length):
		rtn.append(random.uniform(-1,1))
	return rtn
	
def generateRandomSeismogram(length):
	"""
	generates a random 3 compomnent seismogram of size length.
	"""
	ds = []
	for index in range(3):
		ds.append(randomList(length))
	return ds
	
	
def plotSeismograms(dTarget, dOptimized, counter, stationDict):
	"""
	Plot a set of seismograms
	"""
	stationNumber = len(dTarget)
	fig, axes = plt.subplots(stationNumber, 3, sharey=True)
	components = ["z", "r", "t"]
	for vIndex in range(stationNumber):
		for index in range(3):
			axes[vIndex][index].plot(dTarget[vIndex][index])
			axes[vIndex][index].plot(dOptimized[vIndex][index])
			if (vIndex == 0):
				axes[vIndex][index].set_title(components[index])
				
			if (vIndex == (stationNumber) - 1):
				axes[vIndex][index].set_xlabel("t", fontsize=10)
			if (index == 0):
				axes[vIndex][index].set_ylabel(list(stationDict.keys())[vIndex], fontsize=5)
			

	plt.savefig("plots/optimized"+str(counter)+".png")
	
	print("done plotting")
	
	
def main():
	momentTensorScale = 1 ## scale of the moment tensors elements
	popSize = 100 ## size of the genetic population, should be made very large > 1000 for good results
	epochs = 50 ## number of generations the genetic algorithm can go through
	geneticNoise = 1.5 ## scale of random noise within the genetic algorithm. For values less than 2 convergence to some value is garenteed.
	
	
	targetMT = [[0.25, 0, 0], [0, 0.25, 0], [0,0,-0.5]]
	
	##targetMT = randomMT(1)
	depth = 5
	
	for index in range(3):
		print(targetMT[index])
	
	g, stationDict = readGreensFunctions("Greensfunctions", "stationInformation.csv", depth)
	
	## We are not using real data, only testing the programs ability to reproduce 
	## a target moment tensor
	targetData = forwardsProblem(targetMT, g, stationDict) 
	print("done generating target data")
	
	
	optimizedMT = runGA(targetMT, momentTensorScale, popSize, epochs, geneticNoise, g, stationDict)	

	
	
	
main()
	
	
	
	
	
	
	
	

	
	
	

	
	
	
		
		
	
	
	
	


	
