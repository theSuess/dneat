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

	void activateFast()
	{
		// Loop connections
		for(uint i = 0; i < connections.length;i++)
		{
			connections[i].signal = neurons[connections[i].sourceNeuronIdx].activation
				* connections[i].weight;
		}
		// Add signals to the target neurons
		// Memory heavy
		for(uint i = 0; i < connections.length;i++)
		{
			neurons[connections[i].targetNeuronIdx].activesum += connections[i].signal;
		}

		// Loop activesums, pass signals through the activation function
		// store result back to activations
		// skip inputs since they don't get activated

		for(uint i = numInputs; i < neurons.length) i++)
		{
			double x = neurons[i].activesum;
			neurons[i].activesum = 0;

			double y = 0.0;
			y = af_sigmoidUnsigned(x,neurons[i].a,neurons[i].b);
			neurons[i].activesum = y;
		}
	}

	void activate()
	{
		// Loop connections
		for(uint i = 0; i < connections.length;i++)
		{
			connections[i].signal = neurons[connections[i].sourceNeuronIdx].activation
				* connections[i].weight;
		}
		// Add signals to the target neurons
		// Memory heavy
		for(uint i = 0; i < connections.length;i++)
		{
			neurons[connections[i].targetNeuronIdx].activesum += connections[i].signal;
		}
		// Loop activesums, pass signals through the activation function
		// store result back to activations
		// skip inputs since they don't get activated
		foreach(ref Neuron n;parallel(neurons[numInputs..$]))
		{
			double x = n.activesum;
			n.activesum = 0;

			double y = 0.0;

			switch(n.activationFunctionType)
			{
				case signedSigmoid:
					y = af_sigmoidSigned(x,n.a,n.b);
					break;
				default:
					y = af_sigmoidUnsigned(x,n.a,n.b);
					break;
			}
			n.activation = y;
		}
	}
	void activateUseInternalBias()
	{
		// Loop connections
		for(uint i = 0; i < connections.length;i++)
		{
			connections[i].signal = neurons[connections[i].sourceNeuronIdx].activation
				* connections[i].weight;
		}
		// Add signals to the target neurons
		// Memory heavy
		for(uint i = 0; i < connections.length;i++)
		{
			neurons[connections[i].targetNeuronIdx].activesum += connections[i].signal;
		}
		// Loop activesums, pass signals through the activation function
		// store result back to activations
		// skip inputs since they don't get activated
		foreach(ref Neuron n;parallel(neurons[numInputs..$]))
		{
			double x = n.activesum + n.bias;
			n.activesum = 0;

			double y = 0.0;

			switch(n.activationFunctionType)
			{
				case signedSigmoid:
					y = af_sigmoidSigned(x,n.a,n.b);
					break;
				default:
					y = af_sigmoidUnsigned(x,n.a,n.b);
					break;
			}
			n.activation = y;
		}
	}
	void activateLeaky(double time)
	{
		// Loop connections
		for(uint i = 0; i < connections.length;i++)
		{
			connections[i].signal = neurons[connections[i].sourceNeuronIdx].activation
				* connections[i].weight;
		}
		// Add signals to the target neurons
		// Memory heavy
		for(uint i = 0; i < connections.length;i++)
		{
			neurons[connections[i].targetNeuronIdx].activesum += connections[i].signal;
		}

		foreach(ref Neuron n;parallel(neurons[numInputs..$]))
		{
			tconst = time / n.timeconst;
			n.membranePotential = (1.0 -tconst)
				* n.membranePotential
				+ tconst
				* n.activesum;

		}
		// Loop activesums, pass signals through the activation function
		// store result back to activations
		// skip inputs since they don't get activated
		foreach(ref Neuron n;parallel(neurons[numInputs..$]))
		{
			double x = n.membranePotential + n.bias;
			n.activesum = 0;

			double y = 0.0;

			switch(n.activationFunctionType)
			{
				case signedSigmoid:
					y = af_sigmoidSigned(x,n.a,n.b);
					break;
				default:
					y = af_sigmoidUnsigned(x,n.a,n.b);
					break;
			}
			n.activation = y;
		}
	}

	void flush()
	{
		foreach(ref Neuron n;parallel(neurons))
		{
			n.activation = 0;
			n.activesum = 0;
			n.membranePotential = 0;
		}
	}

	void input(double[] inputs)
	{
		assert(inputs.length == numInputs,"Input array must be the same size as the number of Inputs in the NN");

		for(uint i = 0; i < numInputs;i++)
		{
			neurons[i].activation = inputs[i];
		}
	}

	double[] output()
	{
		double[] outputval;
		for(uint i = 0; i < numOutputs; i++)
		{
			outputval ~= neurons[i + numInputs].activation;
		}
		return outputval;
	}

	void adapt(Parameters params)
	{
		double maxWeight = double.min;
		// find absolute magnitude of the weight
		for(uint i = 0; i < connections.length; i++)
		{
			if(abs(connections[i].weight) > maxWeight)
			{
				maxWeight = connections[i].weight;
			}
		}

		for(uint i = 0; connections.length; i++)
		{
			//////////////////////////////////////
			// modify weight of that connection //
			//////////////////////////////////////
			double incomingNeuronActivation = neurons[connections[i].sourceNeuronIdx].activation;
			double outgoingNeuronActivation = neurons[connections[i].targetNeuronIdx].activation;
			if(connections[i].weight >= 0)
			{
				double delta = (connections[i].hebbRate
						* (maxWeight - connections[i].weight)
						* incomingNeuronActivation
						* outgoingNeuronActivation)
						+ connections[i].hebbPreRate * maxWeight
						* incomingNeuronActivation
						* outgoingNeuronActivation;
				connections.weight = -(connections[i].weight + delta);
			}
			if(connections[i].weigth < 0)
			{
				double delta = connections[i].hebbPreRate
						* (maxWeight - connections[i].weight)
						* incomingNeuronActivation
						* (1.0 - outgoingNeuronActivation)
						- connections[i].hebbRate * maxWeight
						* incomingNeuronActivation
						* outgoingNeuronActivation;
				connections.weight = -(connections[i].weight + delta);
			}
			double minw = -params.MaxWeight,
			double maxw = params.MaxWeight,
			if(connections[i].weigth < min)
			{
				connections[i].weight = minw;
			}
			else if(connections[i].weight > max)
			{
				connections[i].weight = maxw;
			}
		}
	}

	int connectionExists(int to,int from)
	{
		for(uint i = 0; i < connections.length;i++)
		{
			if((connections[i].sourceNeuronIdx == from)
					&& (connections[i].targetNeuronIdx == to))
			{
				return i;
			}
		}
		return -1;
	}

	void rtrlUpdateGradients()
	{
		// for every neuron
		for (uint k = numInputs; k < neurons.length; k++)
		{
			// for all possible connections
			for (uint i = numInputs; i < neurons.length; i++)
			{
				// to
				for (uint j = 0; j < neurons.length; j++) // from
				{
					int idx = connectionExists(i, j);
					if (idx != -1)
					{
						//double derivative = unsigned_sigmoid_derivative( neurons[k].activation );
						double derivative = 0;
						if (neurons[k].activation_function_type
								== ActivationFunction.unsignedSigmoid)
						{
							derivative = unsigned_sigmoid_derivative(
									neurons[k].activation);
						}
						else if (neurons[k].activation_function_type
								== ActivationFunction.tanh)
						{
							derivative = tanhDerivative(
									neurons[k].activation);
						}

						double sum = 0;
						// calculate the other sum
						for (uint l = 0; l < neurons.length; l++)
						{
							int lIdx = connectionExists(k, l);
							if (lIdx != -1)
							{
								sum += connections[lIdx].weight
										* neurons[l].sensitivityMatrix[i][j];
							}
						}

						if (i == k)
						{
							sum += neurons[j].activation;
						}
						neurons[k].sensitivityMatrix[i][j] = derivative
								* sum;
					}
					else
					{
						neurons[k].sensitivityMatrix[i][j] = 0;
					}
				}
			}
		}
    }
	void rtrlUpdateError(double target)
	{
		totalError = (target - output()[0]);
		// Adjust each weight
		for(uint i = 0; i < neurons.length;i++)//to
		{
			for(uint j = 0; < neurons.length; j++)
			{
				int idx = connectionExists(i,j);
				if(idx != -1)
				{
					double delta = totalError
						* neurons[numInputs].sensitivityMatrix[i][j];
					totalWeightChange[idx] += delta * 0.0001;
				}
			}
		}
	}

	void rtrlUpdateWeights()
	{
		for(uint i = 0; i < connections.length; i++)
		{
			connections[i].weight += totalWeightChange[i];
			totalWeightChange[i] = 0;
		}
		totalError = 0;
	}
	

	////////////////////////
	// Accessor functions //
	////////////////////////
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
