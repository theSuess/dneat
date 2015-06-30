module dneat.neuralnetwork;

public import dneat.neuralnetwork.connection;
public import dneat.neuralnetwork.neuron;

class NeuralNetwork
{
	double totalError;

	double[] totalWeightChange; // Array size = Number of connections

public:
	ushort numInputs,numOutputs;
	Connection[] connections; // Array size = Number of connections
	Neuron[] neurons;


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
