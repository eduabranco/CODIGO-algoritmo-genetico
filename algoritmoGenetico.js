class GeneticAlgorithm {
    constructor() {
        this.populationSize = 20;
        this.maxGenerations = 1000;
        this.crossoverRate = 0.8;
        this.mutationRate = 0.03;
        this.elitism = true;
        this.boardSize = 8;
        this.individualLength = 24; // 8 rainhas, 3 bits cada (binário)
        this.currentGeneration = 0; // Contador de gerações reais
    }

    // Função para gerar a população inicial
    generatePopulation() {
        let population = [];
        for (let i = 0; i < this.populationSize; i++) {
            population.push(this.generateIndividual());
        }
        return population;
    }

    // Gerar um indivíduo aleatório (binário)
    generateIndividual() {
        let individual = "";
        for (let i = 0; i < this.individualLength; i++) {
            individual += Math.random() < 0.5 ? '0' : '1';
        }
        return individual;
    }

    // Decodificar um indivíduo
    decodeIndividual(individual) {
        let queens = [];
        for (let i = 0; i < this.boardSize; i++) {
            let position = parseInt(individual.slice(i * 3, (i + 1) * 3), 2);
            queens.push(position);
        }
        return queens;
    }

    // Avaliação de fitness (contar conflitos entre rainhas)
    fitness(individual) {
        let queens = this.decodeIndividual(individual);
        let conflicts = 0;

        for (let i = 0; i < queens.length; i++) {
            for (let j = i + 1; j < queens.length; j++) {
                if (queens[i] === queens[j] || Math.abs(queens[i] - queens[j]) === Math.abs(i - j)) {
                    conflicts++;
                }
            }
        }
        return conflicts; // Fitness é o número de conflitos. 0 é o melhor.
    }

    // Seleção por roleta
    rouletteSelection(population) {
        let fitnessSum = population.reduce((sum, individual) => sum + (28 - this.fitness(individual)), 0); // Ajuste para fitness invertido
        let randomValue = Math.random() * fitnessSum;
        let partialSum = 0;

        for (let individual of population) {
            partialSum += (28 - this.fitness(individual));
            if (partialSum >= randomValue) {
                return individual;
            }
        }
    }

    // Cruzamento de um ponto de corte
    crossover(parent1, parent2) {
        if (Math.random() < this.crossoverRate) {
            let point = Math.floor(Math.random() * this.individualLength);
            return [
                parent1.slice(0, point) + parent2.slice(point),
                parent2.slice(0, point) + parent1.slice(point)
            ];
        } else {
            return [parent1, parent2];
        }
    }

    // Mutação (bit flip)
    mutate(individual) {
        let mutated = "";
        for (let i = 0; i < this.individualLength; i++) {
            if (Math.random() < this.mutationRate) {
                mutated += individual[i] === '0' ? '1' : '0';
            } else {
                mutated += individual[i];
            }
        }
        return mutated;
    }

    // Execução do algoritmo genético
    run() {
        let population = this.generatePopulation();
        let bestSolution = null;
        let bestFitness = Infinity;
        this.currentGeneration = 0;

        for (let generation = 0; generation < this.maxGenerations; generation++) {
            this.currentGeneration = generation;
            population = this.evolve(population);

            let currentBest = population.reduce((best, individual) => {
                let fitness = this.fitness(individual);
                if (fitness < bestFitness) {
                    bestSolution = individual;
                    bestFitness = fitness;
                }
                return fitness < this.fitness(best) ? individual : best;
            }, population[0]);

            if (this.fitness(bestSolution) === 0) {
                break; // Encontrou a solução perfeita
            }
        }

        return { bestSolution, bestFitness, generations: this.currentGeneration };
    }

    // Evolução da população
    evolve(population) {
        let newPopulation = [];

        // Elitismo: preserva os melhores indivíduos
        if (this.elitism) {
            newPopulation.push(this.getBest(population));
        }

        // Seleção, cruzamento e mutação
        while (newPopulation.length < this.populationSize) {
            let parent1 = this.rouletteSelection(population);
            let parent2 = this.rouletteSelection(population);
            let [child1, child2] = this.crossover(parent1, parent2);

            newPopulation.push(this.mutate(child1));
            if (newPopulation.length < this.populationSize) {
                newPopulation.push(this.mutate(child2));
            }
        }

        return newPopulation;
    }

    // Pegar o melhor indivíduo
    getBest(population) {
        return population.reduce((best, individual) => this.fitness(individual) < this.fitness(best) ? individual : best, population[0]);
    }
}

// Função para executar o algoritmo 50 vezes e coletar estatísticas
function runExperiments() {
    let iterations = [];
    let times = [];
    let solutions = new Map(); // Usar Map para armazenar solução e fitness

    for (let i = 0; i < 50; i++) {
        let startTime = performance.now();
        let ga = new GeneticAlgorithm();
        let { bestSolution, bestFitness, generations } = ga.run();
        let endTime = performance.now();

        iterations.push(generations); // Número real de gerações
        times.push(endTime - startTime); // Tempo de execução

        solutions.set(bestSolution, bestFitness); // Armazena a solução e o fitness
    }

    // Calcula média e desvio padrão
    let avgIterations = iterations.reduce((sum, val) => sum + val, 0) / iterations.length;
    let avgTime = times.reduce((sum, val) => sum + val, 0) / times.length;
    let stdIterations = Math.sqrt(iterations.reduce((sum, val) => sum + (val - avgIterations) ** 2, 0) / iterations.length);
    let stdTime = Math.sqrt(times.reduce((sum, val) => sum + (val - avgTime) ** 2, 0) / times.length);

    // Exibe os resultados
    console.log(`\nMédia de iterações:\n${avgIterations} \n`);
    console.log(`Desvio padrão de iterações:\n${stdIterations} \n`);
    console.log(`Média de tempo:\n${avgTime} ms \n`);
    console.log(`Desvio padrão de tempo:\n${stdTime} ms \n`);
    console.log(`\nCinco melhores soluções distintas e suas heurísticas (conflitos):\n`);
    
    // Ordenar soluções pelo fitness (número de conflitos) e exibir as 5 melhores
    let sortedSolutions = [...solutions.entries()].sort((a, b) => a[1] - b[1]).slice(0, 5);
    sortedSolutions.forEach(([solution, fitness], index) => {
        console.log(`Solução ${index + 1}: ${solution} | Conflitos: ${fitness}`);
    });
}

// Executar os experimentos
runExperiments();
