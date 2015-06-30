module dneat.genes.neuronGene;

import dneat.genes;

class NeuronGene
{
private:

	// Its unique identification number
	uint _ID;

	// Its type and role in the network
	NeuronType _type;

public:
	// Safe to access directly
	// Modified during evolution

	// useful for displaying the genome
	int x,y;

	// Position (depth) within the network
	double splitY;

	// Additional parameters associated with the
	// neuron's activation function.
	// The current activation function may not use
	// any of them anyway.
	// A is usually used to alter the function's slope with a scalar
	// B is usually used to force a bias to the neuron
	// -------------------
	// Sigmoid : using A, B (slope, shift)
	// Step    : using B    (shift)
	// Gauss   : using A, B (slope, shift))
	// Abs     : using B    (shift)
	// Sine    : using A    (frequency, phase)
	// Square  : using A, B (high phase lenght, low phase length)
	// Linear  : using B    (shift)
	double a,b = 0;

	// Time constant value used when the neuron is activating in leaky integrator mode
	double timeConstant = 0;

	// Bias value used when the neuron is activating in leaky integrator mode
	double bias = 0;
	
	ActivationFunction actFunction = ActivationFunction.signedSigmoid;

	// Constructors
	this(NeuronType type, uint id, double splity)
	{
		_type = type;
		_ID = id;
		this.splitY = splity;
	}

	// Getters for private variables

	uint ID() const @property
	{
		return _ID;
	}

	NeuronType type() const @property
	{
		return _type;
	}

	// Initializing
	void init(double a,double b,double time,double bias,ActivationFunction act)
	{
		this.a = a;
		this.b = b;
		this.timeConstant = time;
		this.bias = bias;
		this.actFunction = act;
	}	
}
