buttons = {
		"A", 
		"B", 
		"Up", 
		"Right", 
		"Down", 
		"Left", 
	}
screenSize = 6
inputSize = (screenSize * 2 + 1) * (screenSize * 2 + 1)

inputs = inputSize + 1
outputs = #buttons

stagnant = 15
increase = 0.1
stopMutations = 0.4
populationX = 300
disjunction = 2.0
weight = 0.4
limit = 1.0

timer = 20

nodes = 1000000

startMutations = 0.2
mutationChance = 0.25
changeChance = 0.9
crossChance = 0.75
linkChance = 2.0
nodeChance = 0.5
biasChance = 0.4

function getMarioPositions()
	positionX = memory.readbyte(0x6D) * 0x100 + memory.readbyte(0x86)
	positionY = memory.readbyte(0x03B8) + 16
	
	screenX = memory.readbyte(0x03AD)
	screenY = memory.readbyte(0x03B8)

end

function getGameSprites()
	local gameSprites = {}
	for sprite = 0, 4 do
		local goomba = memory.readbyte(0xF + sprite)
		if goomba ~=  0 then
			local goombaX = memory.readbyte(0x6E + sprite) * 0x100 + memory.readbyte(0x87 + sprite)
			local goombaY = memory.readbyte(0xCF + sprite)+24
			gameSprites[#gameSprites + 1] = {["x"] = goombaX, ["y"] = goombaY}
		end
	end
	
	return gameSprites
end

function getGameTiles(x, y)
	local gameX = positionX + x + 8
	local gameY = positionY + y - 16
	local gameScreen = math.floor(gameX / 256)%2
	local gameA = math.floor((gameX % 256) / 16)
	local gameB = math.floor((gameY - 32) / 16)
	local addr = 0x500 + gameScreen * 13 * 16 + gameB * 16 + gameA
	
	if gameB >=  13 or gameB < 0 then
		return 0
	end
	
	if memory.readbyte(addr) ~=  0 then
		return 1
	else
		return 0
	end
end
function getinputs()
	getMarioPositions()
	
	gameSprites = getGameSprites()
	screenSprites = {}
	
	local gameinputs = {}
	
	for y =-screenSize * 16, screenSize * 16, 16 do
		for x =-screenSize * 16, screenSize * 16, 16 do
			gameinputs[#gameinputs + 1] = 0
			
			gameTile = getGameTiles(x, y)
			if gameTile ==  1 and positionY+y < 0x1B0 then
				gameinputs[#gameinputs] = 1
			end
			
			for i = 1, #gameSprites do
				distx = math.abs(gameSprites[i]["x"] - (positionX+x))
				disty = math.abs(gameSprites[i]["y"] - (positionY+y))
				if distx <=  8 and disty <=  8 then
					gameinputs[#gameinputs] =  -1
				end
			end

			for i = 1, #screenSprites do
				distx = math.abs(screenSprites[i]["x"] - (positionX+x))
				disty = math.abs(screenSprites[i]["y"] - (positionY+y))
				if distx < 8 and disty < 8 then
					gameinputs[#gameinputs] = -1
				end
			end
		end
	end
	return gameinputs
end

function learningCurve(x)
	return 2 / (1 + math.exp(-4.9 * x)) - 1
end


function newGrouping()
	local grouping = {}
	grouping.species = {}
	grouping.generation = 0
	grouping.group = outputs
	grouping.currentSpecies = 1
	grouping.currentGenome = 1
	grouping.currentFrame = 0
	grouping.bestFitness = 0
	
	return grouping
end

function newGroup()
	grouping.group = grouping.group + 1
	return grouping.group
end

function newSpecies()
	local species = {}
	species.bestFitness = 0
	species.stagnation = 0
	species.genomes = {}
	species.averageFitness = 0
	return species
end

function newGenome()
	local genome = {}
	genome.genes = {}
	genome.fitness = 0
	genome.adjustedFitness = 0
	genome.neuralNetwork = {}
	genome.maxneuron = 0
	genome.globalRank = 0
	genome.mutationRates = {}
	genome.mutationRates["mutationChance"] = mutationChance
	genome.mutationRates["linkChance"] = linkChance
	genome.mutationRates["biasChance"] = biasChance
	genome.mutationRates["nodeChance"] = nodeChance
	genome.mutationRates["start"] = startMutations
	genome.mutationRates["stop"] = stopMutations
	genome.mutationRates["increase"] = increase
	
	return genome
end

function copyGenome(genome)
	local genome2 = newGenome()
	for x = 1, #genome.genes do
		table.insert(genome2.genes, copyGene(genome.genes[x]))
	end
	genome2.maxneuron = genome.maxneuron
	genome2.mutationRates["mutationChance"] = genome.mutationRates["mutationChance"]
	genome2.mutationRates["linkChance"] = genome.mutationRates["linkChance"]
	genome2.mutationRates["biasChance"] = genome.mutationRates["biasChance"]
	genome2.mutationRates["nodeChance"] = genome.mutationRates["nodeChance"]
	genome2.mutationRates["start"] = genome.mutationRates["start"]
	genome2.mutationRates["stop"] = genome.mutationRates["stop"]
	
	return genome2
end

function genome()
	local genome = newGenome()
	local group = 1

	genome.maxneuron = inputs
	mutate(genome)
	
	return genome
end

function newGene()
	local gene = {}
	gene.In = 0
	gene.out = 0
	gene.weight = 0.0
	gene.enabled = true
	gene.group = 0
	
	return gene
end

function copyGene(gene)
	local gene2 = newGene()
	gene2.In = gene.In
	gene2.out = gene.out
	gene2.weight = gene.weight
	gene2.enabled = gene.enabled
	gene2.group = gene.group
	
	return gene2
end

function newNeuron()
	local neuron = {}
	neuron.In = {}
	neuron.value = 0.0
	
	return neuron
end

function generateNeuralNetwork(genome)
	local neuralNetwork = {}
	neuralNetwork.neurons = {}
	
	for i = 1, inputs do
		neuralNetwork.neurons[i] = newNeuron()
	end
	
	for o = 1, outputs do
		neuralNetwork.neurons[nodes + o] = newNeuron()
	end
	
	table.sort(genome.genes, function (a, b)
		return (a.out < b.out)
	end)
	for i = 1, #genome.genes do
		local gene = genome.genes[i]
		if gene.enabled then
			if neuralNetwork.neurons[gene.out] ==  nil then
				neuralNetwork.neurons[gene.out] = newNeuron()
			end
			local neuron = neuralNetwork.neurons[gene.out]
			table.insert(neuron.In, gene)
			if neuralNetwork.neurons[gene.In] ==  nil then
				neuralNetwork.neurons[gene.In] = newNeuron()
			end
		end
	end
	
	genome.neuralNetwork = neuralNetwork
end

function evaluateNetwork(neuralNetwork, gameinputs)
	table.insert(gameinputs, 1)
	if #gameinputs ~=  inputs then
		console.writeline("Inputs are wrong.")
		return {}
	end
	for i = 1, inputs do
		neuralNetwork.neurons[i].value = gameinputs[i]
	end
	
	for _, neuron in pairs(neuralNetwork.neurons) do
		local sum = 0
		for j = 1, #neuron.In do
			local incoming = neuron.In[j]
			local other = neuralNetwork.neurons[incoming.In]
			sum = sum + incoming.weight * other.value
		end
		
		if #neuron.In > 0 then
			neuron.value = learningCurve(sum)
		end
	end
	
	local gameOutputs = {}


	for o = 1, outputs do
		local button = "P1 " .. buttons[o]
		if neuralNetwork.neurons[nodes + o].value > 0 then
			gameOutputs[button] = true
		else
			gameOutputs[button] = false
		end
	end
	
	return gameOutputs
end

function crossover(geneX, geneY)
	if geneY.fitness > geneX.fitness then
		tempGene = geneX
		geneX = geneY
		geneY = tempGene
	end

	local child = newGenome()
	
	local groups2 = {}
	for i = 1, #geneY.genes do
		local gene = geneY.genes[i]
		groups2[gene.group] = gene
	end
	
	for i = 1, #geneX.genes do
		local gene1 = geneX.genes[i]
		local gene2 = groups2[gene1.group]
		if gene2 ~=  nil and math.random(2) ==  1 and gene2.enabled then
			table.insert(child.genes, copyGene(gene2))
		else
			table.insert(child.genes, copyGene(gene1))
		end
	end
	
	child.maxneuron = math.max(geneX.maxneuron, geneY.maxneuron)
	
	for mutation, rate in pairs(geneX.mutationRates) do
		child.mutationRates[mutation] = rate
	end
	
	return child
end

function randomNeuron(genes, nonInput)
	local neurons = {}
	if not nonInput then
		for i = 1, inputs do
			neurons[i] = true
		end
	end
	for o = 1, outputs do
		neurons[nodes+o] = true
	end
	for i = 1, #genes do
		if (not nonInput) or genes[i].In > inputs then
			neurons[genes[i].In] = true
		end
		if (not nonInput) or genes[i].out > inputs then
			neurons[genes[i].out] = true
		end
	end

	local count = 0
	for _, _ in pairs(neurons) do
		count = count + 1
	end
	local n = math.random(1, count)
	
	for k, v in pairs(neurons) do
		n = n - 1
		if n ==  0 then
			return k
		end
	end
	
	return 0
end

function containsLink(genes, link)
	for i = 1, #genes do
		local gene = genes[i]
		if gene.In ==  link.In and gene.out ==  link.out then
			return true
		end
	end
end

function disjoint(genes1, genes2)
	local array1 = {}
	for i = 1, #genes1 do
		local gene = genes1[i]
		array1[gene.group] = true
	end

	local array2 = {}
	for i = 1, #genes2 do
		local gene = genes2[i]
		array2[gene.group] = true
	end
	
	local disjointGenes = 0
	for i = 1, #genes1 do
		local gene = genes1[i]
		if not array2[gene.group] then
			disjointGenes = disjointGenes+1
		end
	end
	
	for i = 1, #genes2 do
		local gene = genes2[i]
		if not array1[gene.group] then
			disjointGenes = disjointGenes+1
		end
	end
	
	local n = math.max(#genes1, #genes2)
	
	return disjointGenes / n
end

function mutate(genome)
	for mutation, rate in pairs(genome.mutationRates) do
		if math.random(1, 2) ==  1 then
			genome.mutationRates[mutation] = 0.95 * rate
		else
			genome.mutationRates[mutation] = 1.05263 * rate
		end
	end
	if math.random() < genome.mutationRates["mutationChance"] then
		pointMutation(genome)
	end
	
	local p = genome.mutationRates["linkChance"]
	while p > 0 do
		if math.random() < p then
			linkMutatione(genome, false)
		end
		p = p - 1
	end

	p = genome.mutationRates["biasChance"]
	while p > 0 do
		if math.random() < p then
			linkMutatione(genome, true)
		end
		p = p - 1
	end
	
	p = genome.mutationRates["nodeChance"]
	while p > 0 do
		if math.random() < p then
			nodeMutation(genome)
		end
		p = p - 1
	end
	
	p = genome.mutationRates["start"]
	while p > 0 do
		if math.random() < p then
			startStopMutation(genome, true)
		end
		p = p - 1
	end

	p = genome.mutationRates["stop"]
	while p > 0 do
		if math.random() < p then
			startStopMutation(genome, false)
		end
		p = p - 1
	end
end

function pointMutation(genome)
	local step = genome.mutationRates["increase"]
	for i = 1, #genome.genes do
		local gene = genome.genes[i]
		if math.random() < changeChance then
			gene.weight = gene.weight + math.random() * increase * 2 - increase
		else
			gene.weight = math.random() * 4 - 2
		end
	end
end

function linkMutatione(genome, forceBias)
	local neuron1 = randomNeuron(genome.genes, false)
	local neuron2 = randomNeuron(genome.genes, true)
	 
	local newLink = newGene()
	if neuron1 <=  inputs and neuron2 <=  inputs then
		return
	end
	if neuron2 <=  inputs then
		local temp = neuron1
		neuron1 = neuron2
		neuron2 = temp
	end

	newLink.In = neuron1
	newLink.out = neuron2
	if forceBias then
		newLink.In = inputs
	end
	
	if containsLink(genome.genes, newLink) then
		return
	end
	newLink.group = newGroup()
	newLink.weight = math.random() * 4 - 2
	
	table.insert(genome.genes, newLink)
end

function nodeMutation(genome)
	if #genome.genes ==  0 then
		return
	end

	genome.maxneuron = genome.maxneuron + 1

	local gene = genome.genes[math.random(1, #genome.genes)]
	if not gene.enabled then
		return
	end
	gene.enabled = false
	
	local gene1 = copyGene(gene)
	gene1.out = genome.maxneuron
	gene1.weight = 1.0
	gene1.group = newGroup()
	gene1.enabled = true
	table.insert(genome.genes, gene1)
	
	local gene2 = copyGene(gene)
	gene2.In = genome.maxneuron
	gene2.group = newGroup()
	gene2.enabled = true
	table.insert(genome.genes, gene2)
end

function startStopMutation(genome, enable)
	local candidates = {}
	for _, gene in pairs(genome.genes) do
		if gene.enabled ==  not enable then
			table.insert(candidates, gene)
		end
	end
	
	if #candidates ==  0 then
		return
	end
	
	local gene = candidates[math.random(1, #candidates)]
	gene.enabled = not gene.enabled
end

function weights(genes1, genes2)
	local array2 = {}
	for i = 1, #genes2 do
		local gene = genes2[i]
		array2[gene.group] = gene
	end

	local sum = 0
	local coincident = 0
	for i = 1, #genes1 do
		local gene = genes1[i]
		if array2[gene.group] ~=  nil then
			local gene2 = array2[gene.group]
			sum = sum + math.abs(gene.weight - gene2.weight)
			coincident = coincident + 1
		end
	end
	
	return sum / coincident
end
	
function sameSpecies(genome1, genome2)
	local disjointMath =  disjunction * disjoint(genome1.genes, genome2.genes)
	local weightMath = weight * weights(genome1.genes, genome2.genes) 
	return disjointMath + weightMath < limit
end

function rankGlobally()
	local global = {}
	for s = 1, #grouping.species do
		local species = grouping.species[s]
		for g = 1, #species.genomes do
			table.insert(global, species.genomes[g])
		end
	end
	table.sort(global, function (a, b)
		return (a.fitness < b.fitness)
	end)
	
	for g = 1, #global do
		global[g].globalRank = g
	end
end

function calculateAverageFitness(species)
	local total = 0
	
	for g = 1, #species.genomes do
		local genome = species.genomes[g]
		total = total + genome.globalRank
	end
	
	species.averageFitness = total / #species.genomes
end

function totalAverageFitness()
	local total = 0
	for s = 1, #grouping.species do
		local species = grouping.species[s]
		total = total + species.averageFitness
	end

	return total
end

function cullSpecies(highLander)
	for s = 1, #grouping.species do
		local species = grouping.species[s]
		
		table.sort(species.genomes, function (a, b)
			return (a.fitness > b.fitness)
		end)
		
		local remaining = math.ceil(#species.genomes / 2)
		if highLander then
			remaining = 1
		end
		while #species.genomes > remaining do
			table.remove(species.genomes)
		end
	end
end

function breedChild(species)
	local child = {}
	if math.random() < crossChance then
		geneX = species.genomes[math.random(1, #species.genomes)]
		geneY = species.genomes[math.random(1, #species.genomes)]
		child = crossover(geneX, geneY)
	else
		g = species.genomes[math.random(1, #species.genomes)]
		child = copyGenome(g)
	end
	
	mutate(child)
	
	return child
end


function removeWeakSpecies()
	local survived = {}

	local sum = totalAverageFitness()
	for s = 1, #grouping.species do
		local species = grouping.species[s]
		breed = math.floor(species.averageFitness / sum * populationX)
		if breed >=  1 then
			table.insert(survived, species)
		end
	end

	grouping.species = survived
end

function removestagnant()
	local survived = {}

	for s = 1, #grouping.species do
		local species = grouping.species[s]
		
		table.sort(species.genomes, function (a, b)
			return (a.fitness > b.fitness)
		end)
		
		if species.genomes[1].fitness > species.bestFitness then
			species.bestFitness = species.genomes[1].fitness
			species.stagnation = 0
		else
			species.stagnation = species.stagnation + 1
		end
		if species.stagnation < stagnant or species.bestFitness >=  grouping.bestFitness then
			table.insert(survived, species)
		end
	end

	grouping.species = survived
end

function addToSpecies(child)
	local foundSpecies = false
	for s = 1, #grouping.species do
		local species = grouping.species[s]
		if not foundSpecies and sameSpecies(child, species.genomes[1]) then
			table.insert(species.genomes, child)
			foundSpecies = true
		end
	end
	
	if not foundSpecies then
		local childSpecies = newSpecies()
		table.insert(childSpecies.genomes, child)
		table.insert(grouping.species, childSpecies)
	end
end

function newGeneration()
	cullSpecies(false) 
	rankGlobally()
	removestagnant()
	rankGlobally()
	for s = 1, #grouping.species do
		local species = grouping.species[s]
		calculateAverageFitness(species)
	end
	removeWeakSpecies()
	local sum = totalAverageFitness()
	local children = {}
	for s = 1, #grouping.species do
		local species = grouping.species[s]
		breed = math.floor(species.averageFitness / sum * populationX) - 1
		for i = 1, breed do
			table.insert(children, breedChild(species))
		end
	end
	cullSpecies(true)
	while #children + #grouping.species < populationX do
		local species = grouping.species[math.random(1, #grouping.species)]
		table.insert(children, breedChild(species))
	end
	for c = 1, #children do
		local child = children[c]
		addToSpecies(child)
	end
	
	grouping.generation = grouping.generation + 1
	
	
end
	
function startTesting()
	grouping = newGrouping()
	for i = 1, populationX do
		basic = genome()
		addToSpecies(basic)
	end
	startRun()
end

function cleanController()
	controller = {}
	for b = 1, #buttons do
		controller["P1 " .. buttons[b]] = false
	end
	joypad.set(controller)
end

function startRun()
	rightmost = 0
	grouping.currentFrame = 0
	timeout = timer
	cleanController()	
	local species = grouping.species[grouping.currentSpecies]
	local genome = species.genomes[grouping.currentGenome]
	generateNeuralNetwork(genome)
	testCurrent()
end

function testCurrent()
	local species = grouping.species[grouping.currentSpecies]
	local genome = species.genomes[grouping.currentGenome]

	gameinputs = getinputs()
	controller = evaluateNetwork(genome.neuralNetwork, gameinputs)
	if controller["P1 Left"] and controller["P1 Right"] then
		controller["P1 Left"] = false
		controller["P1 Right"] = false
	end
	if controller["P1 Up"] and controller["P1 Down"] then
		controller["P1 Up"] = false
		controller["P1 Down"] = false
	end

	joypad.set(controller)
end

if grouping ==  nil then
	startTesting()
end

function nextGenome()
	grouping.currentGenome = grouping.currentGenome + 1
	if grouping.currentGenome > #grouping.species[grouping.currentSpecies].genomes then
		grouping.currentGenome = 1
		grouping.currentSpecies = grouping.currentSpecies+1
		if grouping.currentSpecies > #grouping.species then
			newGeneration()
			grouping.currentSpecies = 1
		end
	end
end

function gotEm()
	local species = grouping.species[grouping.currentSpecies]
	local genome = species.genomes[grouping.currentGenome]
	
	return genome.fitness ~=  0
end

function showButtons(genome)
	local neuralNetwork = genome.neuralNetwork
	local cells = {}
	local cell = {}
	for o = 1, outputs do
		cell = {}
		cell.gameX = 220
		cell.gameY = 30 + 8 * o
		cell.value = neuralNetwork.neurons[nodes + o].value
		cells[nodes+o] = cell
		local color
		if cell.value > 0 then
			color = 0xFF0000FF
		else
			color = 0xFF000000
		end
		gui.drawText(0, 140+8 * o, buttons[o], color, 9)
	end
end

while gotEm() do
		extGenome()
end
startRun()
grouping.currentFrame = grouping.currentFrame + 1

while true do
	memory.writebyte(0x07FA, 1)
    memory.writebyte(0x075A, 3)
	local species = grouping.species[grouping.currentSpecies]

	local genome = species.genomes[grouping.currentGenome]
	
	showButtons(genome)
	
	if grouping.currentFrame % 5 ==  0 then
		testCurrent()
	end

	joypad.set(controller)

	getMarioPositions()
	if positionX > rightmost then
		rightmost = positionX
		timeout = timer
	end
	
	timeout = timeout - 1
	
	
	local timeoutBonus = grouping.currentFrame / 4
	if timeout + timeoutBonus <=  0 then
		local fitness = rightmost - grouping.currentFrame / 2
		if rightmost > 3186 then
			fitness = fitness + 1000
		end
		if fitness ==  0 then
			fitness =  -1
		end
		genome.fitness = fitness
		
		if fitness > grouping.bestFitness then
			grouping.bestFitness = fitness
		end
		
		grouping.currentSpecies = 1
		grouping.currentGenome = 1
		while gotEm() do
			nextGenome()
		end
		startRun()
	end

	local measured = 0
	local total = 0

	for _, species in pairs(grouping.species) do
		for _, genome in pairs(species.genomes) do
			total = total + 1
			if genome.fitness ~=  0 then
				measured = measured + 1
			end
		end
	end
		gui.drawText(0, 36, "Distance: " .. math.floor(rightmost - (grouping.currentFrame) / 2 - (timeout + timeoutBonus) * 2 / 3), 0xFF000000, 11)
		gui.drawText(100, 36, "  Best: " .. math.floor(grouping.bestFitness), 0xFF000000, 11)
		gui.drawText(0, 48, "Gen: " .. grouping.generation .. " Species: " .. grouping.currentSpecies .. " Tested: (" .. math.floor(measured / total * 100) .. "%)", 0xFF000000, 30)
	grouping.currentFrame = grouping.currentFrame + 1

	emu.frameadvance();
end
