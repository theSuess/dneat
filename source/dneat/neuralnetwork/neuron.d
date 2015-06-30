module dneat.neuralnetwork.neuron;

import dneat.genes;

class Neuron
{
public:
	// Synaptic input
	double activesum;
	// Synaptic input passed through the activation function
	double activation;

	double a,b,timeconst,bias; // misc params
	double membranePotential; // used in leaky integrator mode

	ActivationFunction activationFunctionType;

	//used for displaying
	double x,y,z;
	double sx,sy,sz;

	double[] substrateCoords;
	double splitY;
	NeuronType type;

	// Used for RTRL learning
	double[][] sensitivityMatrix;

	override bool opEquals(Object o)
	{
		Neuron n = cast(Neuron)o;
		return (type == n.type && splitY == n.splitY && activationFunctionType == n.activationFunctionType);
	}
}
