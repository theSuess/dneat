module dneat.parameters;

//////////////////////////////////////////////
// The NEAT Parameters class
//////////////////////////////////////////////
class Parameters
{
public:
    /////////////////////
    // Members
    /////////////////////


    ////////////////////
    // Basic parameters
    ////////////////////

    // Size of population
    uint PopulationSize;

    // If true, this enables dynamic compatibility thresholding
    // It will keep the number of species between MinSpecies and MaxSpecies
    bool DynamicCompatibility;

    // Minimum number of species
    uint MinSpecies;

    // Maximum number of species
    uint MaxSpecies;

    // Don't wipe the innovation database each generation?
    bool InnovationsForever;

    // Allow clones or nearly identical genomes to exist simultaneously in the population.
    // This is useful for non-deterministic environments,
    // as the same individual will get more than one chance to prove himself, also
    // there will be more chances the same individual to mutate in different ways.
    // The drawback is greatly increased time for reproduction. If you want to
    // search quickly, yet less efficient, leave this to true.
    bool AllowClones;

   ////////////////////////////////
    // GA Parameters
    ////////////////////////////////

    // Age treshold, meaning if a species is below it, it is considered young
    uint YoungAgeTreshold;

    // Fitness boost multiplier for young species (1.0 means no boost)
    // Make sure it is >= 1.0 to avoid confusion
    double YoungAgeFitnessBoost;

    // Number of generations without improvement (stagnation) allowed for a species
    uint SpeciesMaxStagnation;

    // Minimum jump in fitness necessary to be considered as improvement.
    // Setting this value to 0.0 makes the system to behave like regular NEAT.
    double StagnationDelta;

    // Age threshold, meaning if a species if above it, it is considered old
    uint OldAgeTreshold;

    // Multiplier that penalizes old species.
    // Make sure it is < 1.0 to avoid confusion.
    double OldAgePenalty;

    // Detect competetive coevolution stagnation
    // This kills the worst species of age >N (each X generations)
    bool DetectCompetetiveCoevolutionStagnation;

    // Each X generation..
    int KillWorstSpeciesEach;

    // Of age above..
    int KillWorstAge;

    // Percent of best individuals that are allowed to reproduce. 1.0 = 100%
    double SurvivalRate;

    // Probability for a baby to result from sexual reproduction (crossover/mating). 1.0 = 100%
    double CrossoverRate;

    // If a baby results from sexual reproduction, this probability determines if mutation will
    // be performed after crossover. 1.0 = 100% (always mutate after crossover)
    double OverallMutationRate;

    // Probability for a baby to result from inter-species mating.
    double InterspeciesCrossoverRate;

    // Probability for a baby to result from Multipoint Crossover when mating. 1.0 = 100%
    // The default if the Average mating.
    double MultipointCrossoverRate;

    // Performing roulette wheel selection or not?
    bool RouletteWheelSelection;

    ///////////////////////////////////
    // Phased Search parameters   //
    ///////////////////////////////////

    // Using phased search or not
    bool PhasedSearching;

    // Using delta coding or not
    bool DeltaCoding;

    // What is the MPC + base MPC needed to begin simplifying phase
    uint SimplifyingPhaseMPCTreshold;

    // How many generations of global stagnation should have passed to enter simplifying phase
    uint SimplifyingPhaseStagnationTreshold;

    // How many generations of MPC stagnation are needed to turn back on complexifying
    uint ComplexityFloorGenerations;


    /////////////////////////////////////
    // Novelty Search parameters       //
    /////////////////////////////////////

    // the K constant
    uint NoveltySearch_K;

    // Sparseness treshold. Add to the archive if above
    double NoveltySearch_P_min;

    // Dynamic Pmin?
    bool NoveltySearch_Dynamic_Pmin;

    // How many evaluations should pass without adding to the archive
    // in order to lower Pmin
    uint NoveltySearch_No_Archiving_Stagnation_Treshold;

    // How should it be multiplied (make it less than 1.0)
    double NoveltySearch_Pmin_lowering_multiplier;

    // Not lower than this value
    double NoveltySearch_Pmin_min;


    // How many one-after-another additions to the archive should
    // pass in order to raise Pmin
    uint NoveltySearch_Quick_Archiving_Min_Evaluations;

    // How should it be multiplied (make it more than 1.0)
    double NoveltySearch_Pmin_raising_multiplier;

    // Per how many evaluations to recompute the sparseness
    uint NoveltySearch_Recompute_Sparseness_Each;


    ///////////////////////////////////
    // Mutation parameters
    ///////////////////////////////////

    // Probability for a baby to be mutated with the Add-Neuron mutation.
    double MutateAddNeuronProb;

    // Allow splitting of any recurrent links
    bool SplitRecurrent;

    // Allow splitting of looped recurrent links
    bool SplitLoopedRecurrent;

    // Maximum number of tries to find a link to split
    int NeuronTries;

    // Probability for a baby to be mutated with the Add-Link mutation
    double MutateAddLinkProb;

    // Probability for a new incoming link to be from the bias neuron;
    double MutateAddLinkFromBiasProb;

    // Probability for a baby to be mutated with the Remove-Link mutation
    double MutateRemLinkProb;

    // Probability for a baby that a simple neuron will be replaced with a link
    double MutateRemSimpleNeuronProb;

    // Maximum number of tries to find 2 neurons to add/remove a link
    uint LinkTries;

    // Probability that a link mutation will be made recurrent
    double RecurrentProb;

    // Probability that a recurrent link mutation will be looped
    double RecurrentLoopProb;

    // Probability for a baby's weights to be mutated
    double MutateWeightsProb;

    // Probability for a severe (shaking) weight mutation
    double MutateWeightsSevereProb;

    // Probability for a particular gene to be mutated. 1.0 = 100%
    double WeightMutationRate;

    // Maximum perturbation for a weight mutation
    double WeightMutationMaxPower;

    // Maximum magnitude of a replaced weight
    double WeightReplacementMaxPower;

    // Maximum absolute magnitude of a weight
    double MaxWeight;

    // Probability for a baby's A activation function parameters to be perturbed
    double MutateActivationAProb;

    // Probability for a baby's B activation function parameters to be perturbed
    double MutateActivationBProb;

    // Maximum magnitude for the A parameter perturbation
    double ActivationAMutationMaxPower;

    // Maximum magnitude for the B parameter perturbation
    double ActivationBMutationMaxPower;

    // Maximum magnitude for time costants perturbation
    double TimeConstantMutationMaxPower;

    // Maximum magnitude for biases perturbation
    double BiasMutationMaxPower;

    // Activation parameter A min/max
    double MinActivationA;
    double MaxActivationA;

    // Activation parameter B min/max
    double MinActivationB;
    double MaxActivationB;

    // Probability for a baby that an activation function type will be changed for a single neuron
    // considered a structural mutation because of the large impact on fitness
    double MutateNeuronActivationTypeProb;

    // Probabilities for a particular activation function appearance
    double ActivationFunction_SignedSigmoid_Prob;
    double ActivationFunction_UnsignedSigmoid_Prob;
    double ActivationFunction_Tanh_Prob;
    double ActivationFunction_TanhCubic_Prob;
    double ActivationFunction_SignedStep_Prob;
    double ActivationFunction_UnsignedStep_Prob;
    double ActivationFunction_SignedGauss_Prob;
    double ActivationFunction_UnsignedGauss_Prob;
    double ActivationFunction_Abs_Prob;
    double ActivationFunction_SignedSine_Prob;
    double ActivationFunction_UnsignedSine_Prob;
    double ActivationFunction_SignedSquare_Prob;
    double ActivationFunction_UnsignedSquare_Prob;
    double ActivationFunction_Linear_Prob;

    // Probability for a baby's neuron time constant values to be mutated
    double MutateNeuronTimeConstantsProb;

    // Probability for a baby's neuron bias values to be mutated
    double MutateNeuronBiasesProb;

    // Time constant range
    double MinNeuronTimeConstant;
    double MaxNeuronTimeConstant;

    // Bias range
    double MinNeuronBias;
    double MaxNeuronBias;

    /////////////////////////////////////
    // Speciation parameters
    /////////////////////////////////////

    // Percent of disjoint genes importance
    double DisjointCoeff;

    // Percent of excess genes importance
    double ExcessCoeff;

    // Node-specific activation parameter A difference importance
    double ActivationADiffCoeff;

    // Node-specific activation parameter B difference importance
    double ActivationBDiffCoeff;

    // Average weight difference importance
    double WeightDiffCoeff;

    // Average time constant difference importance
    double TimeConstantDiffCoeff;

    // Average bias difference importance
    double BiasDiffCoeff;

    // Activation function type difference importance
    double ActivationFunctionDiffCoeff;

    // Compatibility treshold
    double CompatTreshold;

    // Minumal value of the compatibility treshold
    double MinCompatTreshold;

    // Modifier per generation for keeping the species stable
    double CompatTresholdModifier;

    // Per how many generations to change the treshold
    uint CompatTreshChangeInterval_Generations;

    // Per how many evaluations to change the treshold
    uint CompatTreshChangeInterval_Evaluations;


    // ES HyperNEAT params


    double DivisionThreshold;

    double VarianceThreshold;

    // Used for Band prunning.
    double BandThreshold;

    // Max and Min Depths of the quadtree
    uint InitialDepth;

    uint MaxDepth;

    // How many hidden layers before connecting nodes to output. At 0 there is
    // one hidden layer. At 1, there are two and so on.
    uint IterationLevel;

    // The Bias value for the CPPN queries.
    double CPPN_Bias;

    // Quadtree Dimensions
    // The range of the tree. Typically set to 2,
    double Width;

    // The (x, y) coordinates of the tree
    double Qtree_X;

    double Qtree_Y;

    // Use Link Expression output
    bool Leo;

    // Threshold above which a connection is expressed
    double LeoThreshold;

    // Use geometric seeding. Currently only along the X axis. 1
    bool LeoSeed;

    ////////////////////////////////////
    // Methods
    ////////////////////////////////////
//TODO: Implement
    // Load the parameters from a file
    // returns 0 on success
//int Load(const char* filename);
    // Load the parameters from an already opened file for reading
//int Load(std::ifstream& a_DataFile);

//void Save(const char* filename);
    // Saves the parameters to an already opened file for writing
//void Save(FILE* a_fstream);


    // resets the parameters to built-in defaults
    void Reset()
		{
		////////////////////
		// Basic parameters
		////////////////////

		// Size of population
		PopulationSize = 300;

		// If true, this enables dynamic compatibility thresholding
		// It will keep the number of species between MinSpecies and MaxSpecies
		DynamicCompatibility = true;

		// Minimum number of species
		MinSpecies = 5;

		// Maximum number of species
		MaxSpecies = 10;

		// Don't wipe the innovation database each generation?
		InnovationsForever = true;

		// Allow clones or nearly identical genomes to exist simultaneously in the population.
		// This is useful for non-deterministic environments,
		// as the same individual will get more than one chance to prove himself, also
		// there will be more chances the same individual to mutate in different ways.
		// The drawback is greatly increased time for reproduction. If you want to
		// search quickly, yet less efficient, leave this to true.
		AllowClones = true;




		////////////////////////////////
		// GA Parameters
		////////////////////////////////

		// Age treshold, meaning if a species is below it, it is considered young
		YoungAgeTreshold = 5;

		// Fitness boost multiplier for young species (1.0 means no boost)
		// Make sure it is >= 1.0 to avoid confusion
		YoungAgeFitnessBoost = 1.1;

		// Number of generations without improvement (stagnation) allowed for a species
		SpeciesMaxStagnation = 50;

		// Minimum jump in fitness necessary to be considered as improvement.
		// Setting this value to 0.0 makes the system to behave like regular NEAT.
		StagnationDelta = 0.0;

		// Age threshold, meaning if a species is above it, it is considered old
		OldAgeTreshold = 30;

		// Multiplier that penalizes old species.
		// Make sure it is <= 1.0 to avoid confusion.
		OldAgePenalty = 1.0;

		// Detect competetive coevolution stagnation
		// This kills the worst species of age >N (each X generations)
		DetectCompetetiveCoevolutionStagnation = false;
		// Each X generation..
		KillWorstSpeciesEach = 15;
		// Of age above..
		KillWorstAge = 10;

		// Percent of best individuals that are allowed to reproduce. 1.0 = 100%
		SurvivalRate = 0.2;

		// Probability for a baby to result from sexual reproduction (crossover/mating). 1.0 = 100%
		// If asexual reprodiction is chosen, the baby will be mutated 100%
		CrossoverRate = 0.7;

		// If a baby results from sexual reproduction, this probability determines if mutation will
		// be performed after crossover. 1.0 = 100% (always mutate after crossover)
		OverallMutationRate = 0.25;

		// Probability for a baby to result from inter-species mating.
		InterspeciesCrossoverRate = 0.0001;

		// Probability for a baby to result from Multipoint Crossover when mating. 1.0 = 100%
		// The default is the Average mating.
		MultipointCrossoverRate = 0.75;

		// Performing roulette wheel selection or not?
		RouletteWheelSelection = false;





		///////////////////////////////////
		// Phased Search parameters   //
		///////////////////////////////////

		// Using phased search or not
		PhasedSearching = false;

		// Using delta coding or not
		DeltaCoding = false;

		// What is the MPC + base MPC needed to begin simplifying phase
		SimplifyingPhaseMPCTreshold = 20;

		// How many generations of global stagnation should have passed to enter simplifying phase
		SimplifyingPhaseStagnationTreshold = 30;

		// How many generations of MPC stagnation are needed to turn back on complexifying
		ComplexityFloorGenerations = 40;






		/////////////////////////////////////
		// Novelty Search parameters       //
		/////////////////////////////////////

		// the K constant
		NoveltySearch_K = 15;

		// Sparseness treshold. Add to the archive if above
		NoveltySearch_P_min = 0.5;

		// Dynamic Pmin?
		NoveltySearch_Dynamic_Pmin = true;

		// How many evaluations should pass without adding to the archive
		// in order to lower Pmin
		NoveltySearch_No_Archiving_Stagnation_Treshold = 150;

		// How should it be multiplied (make it less than 1.0)
		NoveltySearch_Pmin_lowering_multiplier = 0.9;

		// Not lower than this value
		NoveltySearch_Pmin_min = 0.05;

		// How many one-after-another additions to the archive should
		// pass in order to raise Pmin
		NoveltySearch_Quick_Archiving_Min_Evaluations = 8;

		// How should it be multiplied (make it more than 1.0)
		NoveltySearch_Pmin_raising_multiplier = 1.1;

		// Per how many evaluations to recompute the sparseness of the population
		NoveltySearch_Recompute_Sparseness_Each = 25;




		///////////////////////////////////
		// Structural Mutation parameters
		///////////////////////////////////

		// Probability for a baby to be mutated with the Add-Neuron mutation.
		MutateAddNeuronProb = 0.01;

		// Allow splitting of any recurrent links
		SplitRecurrent = true;

		// Allow splitting of looped recurrent links
		SplitLoopedRecurrent = true;

		// Probability for a baby to be mutated with the Add-Link mutation
		MutateAddLinkProb = 0.03;

		// Probability for a new incoming link to be from the bias neuron;
		// This enforces it. A value of 0.0 doesn't mean there will not be such links
		MutateAddLinkFromBiasProb = 0.0;

		// Probability for a baby to be mutated with the Remove-Link mutation
		MutateRemLinkProb = 0.0;

		// Probability for a baby that a simple neuron will be replaced with a link
		MutateRemSimpleNeuronProb = 0.0;

		// Maximum number of tries to find 2 neurons to add/remove a link
		LinkTries = 32;

		// Probability that a link mutation will be made recurrent
		RecurrentProb = 0.25;

		// Probability that a recurrent link mutation will be looped
		RecurrentLoopProb = 0.25;





		///////////////////////////////////
		// Parameter Mutation parameters
		///////////////////////////////////

		// Probability for a baby's weights to be mutated
		MutateWeightsProb = 0.90;

		// Probability for a severe (shaking) weight mutation
		MutateWeightsSevereProb = 0.25;

		// Probability for a particular gene's weight to be mutated. 1.0 = 100%
		WeightMutationRate = 1.0;

		// Maximum perturbation for a weight mutation
		WeightMutationMaxPower = 1.0;

		// Maximum magnitude of a replaced weight
		WeightReplacementMaxPower = 1.0;

		// Maximum absolute magnitude of a weight
		MaxWeight = 8.0;

		// Probability for a baby's A activation function parameters to be perturbed
		MutateActivationAProb = 0.0;

		// Probability for a baby's B activation function parameters to be perturbed
		MutateActivationBProb = 0.0;

		// Maximum magnitude for the A parameter perturbation
		ActivationAMutationMaxPower = 0.0;

		// Maximum magnitude for the B parameter perturbation
		ActivationBMutationMaxPower = 0.0;

		// Activation parameter A min/max
		MinActivationA = 1.0;
		MaxActivationA = 1.0;

		// Activation parameter B min/max
		MinActivationB = 0.0;
		MaxActivationB = 0.0;

		// Maximum magnitude for time costants perturbation
		TimeConstantMutationMaxPower = 0.0;

		// Maximum magnitude for biases perturbation
		BiasMutationMaxPower = WeightMutationMaxPower;

		// Probability for a baby's neuron time constant values to be mutated
		MutateNeuronTimeConstantsProb = 0.0;

		// Probability for a baby's neuron bias values to be mutated
		MutateNeuronBiasesProb = 0.0;

		// Time constant range
		MinNeuronTimeConstant = 0.0;
		MaxNeuronTimeConstant = 0.0;

		// Bias range
		MinNeuronBias = 0.0;
		MaxNeuronBias = 0.0;





		// Probability for a baby that an activation function type will be changed for a single neuron
		// considered a structural mutation because of the large impact on fitness
		MutateNeuronActivationTypeProb = 0.0;

		// Probabilities for a particular activation function appearance
		ActivationFunction_SignedSigmoid_Prob = 0.0;
		ActivationFunction_UnsignedSigmoid_Prob = 1.0;
		ActivationFunction_Tanh_Prob = 0.0;
		ActivationFunction_TanhCubic_Prob = 0.0;
		ActivationFunction_SignedStep_Prob = 0.0;
		ActivationFunction_UnsignedStep_Prob = 0.0;
		ActivationFunction_SignedGauss_Prob = 0.0;
		ActivationFunction_UnsignedGauss_Prob = 0.0;
		ActivationFunction_Abs_Prob = 0.0;
		ActivationFunction_SignedSine_Prob = 0.0;
		ActivationFunction_UnsignedSine_Prob = 0.0;
		ActivationFunction_SignedSquare_Prob = 0.0;
		ActivationFunction_UnsignedSquare_Prob = 0.0;
		ActivationFunction_Linear_Prob = 0.0;




		/////////////////////////////////////
		// Speciation parameters
		/////////////////////////////////////

		// Percent of disjoint genes importance
		DisjointCoeff = 1.0;

		// Percent of excess genes importance
		ExcessCoeff = 1.0;

		// Average weight difference importance
		WeightDiffCoeff = 0.5;

		// Node-specific activation parameter A difference importance
		ActivationADiffCoeff = 0.0;

		// Node-specific activation parameter B difference importance
		ActivationBDiffCoeff = 0.0;

		// Average time constant difference importance
		TimeConstantDiffCoeff = 0.0;

		// Average bias difference importance
		BiasDiffCoeff = 0.0;

		// Activation function type difference importance
		ActivationFunctionDiffCoeff = 0.0;

		// Compatibility treshold
		CompatTreshold = 5.0;

		// Minumal value of the compatibility treshold
		MinCompatTreshold = 0.2;

		// Modifier per generation for keeping the species stable
		CompatTresholdModifier = 0.3;

		// Per how many generations to change the treshold
		// (used in generational mode)
		CompatTreshChangeInterval_Generations = 1;

		// Per how many evaluations to change the treshold
		// (used in steady state mode)
		CompatTreshChangeInterval_Evaluations = 10;


		DivisionThreshold = 0.03;

		VarianceThreshold = 0.03;

		// Used for Band prunning.
		BandThreshold = 0.3;

		// Max and Min Depths of the quadtree
		InitialDepth = 3;

		MaxDepth = 3;

		// How many hidden layers before connecting nodes to output. At 0 there is
		// one hidden layer. At 1, there are two and so on.
		IterationLevel = 1;

		// The Bias value for the CPPN queries.
		CPPN_Bias = -1;

		// Quadtree Dimensions
		// The range of the tree. Typically set to 2,
		Width = 2.0;

		// The (x, y) coordinates of the tree
		Qtree_X = 0.0;

		Qtree_Y = 0.0;

		// Use Link Expression output
		Leo = false;

		// Threshold above which a connection is expressed
		LeoThreshold = 0.1;

		// Use geometric seeding. Currently only along the X axis. 1
		LeoSeed = false;	
	}
	this()
	{
		Reset();
	}
}	
