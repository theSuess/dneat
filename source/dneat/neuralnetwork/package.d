module dneat.neuralnetwork;

public import dneat.neuralnetwork.connection;
public import dneat.neuralnetwork.neuron;
import std.math;
import std.parallelism;
import std.random;

class NeuralNetwork
{
	double totalError;

	double[] totalWeightChange; // Array size = Number of connections

public:
	ushort numInputs,numOutputs;
	Connection[] connections; // Array size = Number of connections
	Neuron[] neurons;

	// Constructors
	this(bool minimal)
	{
		if(!minimal)
		{
			//Build an XOR Network

			// Initialize three input neurons
			Neuron i1,i2,i3;

			// One output neuron
			Neuron o;

			// One hidden neuron
			Neuron h;

			// Add the Neurons to the net

			neurons ~= i1;
			neurons ~= i2;
			neurons ~= i3;
			neurons ~= o;
			neurons ~= h;

			// Add the neuron connection

			Connection c = new Connection();
			c.sourceNeuronIdx = 0;
			c.targetNeuronIdx = 3;
			c.weight = 0;
			connections ~= c;

			c.sourceNeuronIdx = 1;
			c.targetNeuronIdx = 3;
			c.weight = 0;
			connections ~= c;

			c.sourceNeuronIdx = 2;
			c.targetNeuronIdx = 3;
			c.weight = 0;
			connections ~= c;

			c.sourceNeuronIdx = 0;
			c.targetNeuronIdx = 4;
			c.weight = 0;
			connections ~= c;

			c.sourceNeuronIdx = 1;
			c.targetNeuronIdx = 4;
			c.weight = 0;
			connections ~= c;

			c.sourceNeuronIdx = 2;
			c.targetNeuronIdx = 4;
			c.weight = 0;
			connections ~= c;

			c.sourceNeuronIdx = 4;
			c.targetNeuronIdx = 3;
			c.weight = 0;
			connections ~= c;

			// Initialize Weights (random)
			foreach(ref Connection conn;parallel(connections))
			{
				conn.weight = uniform(0.0,1)-0.5;
			}

			foreach(ref Neuron n; parallel(neurons))
			{
				n.a = 1;
				n.b = 0;
				n.timeconst = 0;
				n.bias = 0;
				n.membranePotential = 0;
			}
			initRTRLMatrix();
		}
		else
		{
			initEmpty();
		}
	}

	this()
	{
		initEmpty();
	}

	// Init an empty Neural Network
	void initEmpty()
	{
		numInputs = 0;
		numOutputs = 0;
		totalError = 0;
		clear();
	}


	void initRTRLMatrix()
	{
		// Allocate Memory for the sensitivity Matrix and fill it with 0
		foreach(Neuron n; neurons)
		{
			n.sensitivityMatrix.length = neurons.length;
			foreach(ref double[] d; n.sensitivityMatrix)
			{
				d.length = neurons.length;
				foreach(ref double j; d)
				{
					j = 0;
				}
			}
		}
		totalError = 0;
		totalWeightChange = [];
		totalWeightChange.length = connections.length;
		foreach(ref double d; totalWeightChange)
		{
			d = 0;
		}

	}

	// Accessor functions
	void addNeuron(Neuron n)
	{
		neurons ~= n;
	}

	void addConnection(Connection c)
	{
		connections ~= c;
	}

	Connection getConnectionByIndex(uint idx)
	{
		return connections[idx];
	}

	Neuron getNeuronByIndex(uint idx)
	{
		return neurons[idx];
	}

	void setInputOutputDimensions(const ushort i,const ushort o)
	{
		numInputs = i;
		numOutputs = o;
	}

	// Clears the Network and makes it a minimal one
	void clear()
	{
		connections = [];
		neurons = [];
		totalWeightChange = [];
		setInputOutputDimensions(0,0);
	}
}

//////////////////////////
// Activation Functions //
//////////////////////////

double af_sigmoidUnsigned(double x,double slope,double shift)
{
	return 1.0 / (1.0 + exp( - slope * x - shift));
}

double af_sigmoidSigned(double x,double slope,double shift)
{
	double tY = af_sigmoidUnsigned(x,slope,shift);
	return (tY - 0.5) * 2.0;
}

double af_tanh(double x,double slope,double shift)
{
	return tanh(x * slope);
}

double af_tanhCubic(double x,double slope,double shift)
{
	return tanh(x * x * x * slope);
}

double af_stepSigned(double x, double shift)
{
	double tY;
	if (x > shift)
	{
		tY = 1.0;
	}
	else
	{
		tY = -1.0;
	}

	return tY;
}

double af_stepUnsigned(double x, double shift)
{
	if (x > shift)
	{
		return 1.0;
	}
	else
	{
		return 0.0;
	}
}

double af_gaussSigned(double x, double slope, double shift)
{
	double tY = exp( - slope * x * x + shift);
	return (tY-0.5)*2.0;
}

double af_gaussUnsigned(double x, double slope, double shift)
{
    return exp( - slope * x * x + shift);
}

double af_abs(double x, double shift)
{
    return ((x + shift)< 0.0)? -(x + shift): (x + shift);
}

double af_sineSigned(double x, double freq, double shift)
{
	freq = 3.141592;
	return sin(x * freq + shift);
}

double af_sineUnsigned(double x, double freq, double shift)
{
	double tY = sin((x * freq + shift) );
	return (tY + 1.0) / 2.0;
}

double af_squareSigned(double x, double highPulseSize, double lowPulseSize)
{
	return 0.0;    // TODO
}
double af_squareUnsigned(double x, double highPulseSize, double lowPulseSize)
{
	return 0.0;    // TODO
}

double af_linear(double x, double shift)
{
	return (x + shift);
}

double unsignedSigmoidDerivative(double x)
{
	return x * (1 - x);
}

double tanhDerivative(double x)
{
	return 1 - x * x;
}
