module dneat.genes;

public import dneat.genes.linkGene;
public import dneat.genes.neuronGene;

//////////////////////////////////////////////
// Enumeration for all available neuron types
//////////////////////////////////////////////
enum NeuronType
{
	none = 0,
	input,
	bias,
	hidden,
	output
}


//////////////////////////////////////////////////////////
// Enumeration for all possible activation function types
//////////////////////////////////////////////////////////
enum ActivationFunction
{
	signedSigmoid = 0,   // Sigmoid function   (default) (blurred cutting plane)
	unsignedSigmoid,
	tanh,
	tanhCubic,
	signedStep,          // Treshold (0 or 1)  (cutting plane)
	unsignedStep,
	signedGauss,         // Gaussian           (symettry)
	unsignedGauss,
	abs,                  // Absolute value |x| (another symettry)
	signedSine,          // Sine wave          (smooth repetition)
	unsignedSine,
	signedSquare,        // Square wave        (pulse repetition)
	unsignedSquare,
	linear               // Linear f(x)=x      (combining coordinate frames only)
}
