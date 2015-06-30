module dneat.neuralnetwork.connection;

class Connection
{
public:
	// Source and target neuron index
	uint sourceNeuronIdx;
	uint targetNeuronIdx;

	double weight; // weight of the connection
	double signal; // weight * input signal

	// Hebbian learning parameters
	// Ignored in case there is no lifetime learning
	double hebb_rate;
	double hebb_pre_rate;

	override bool opEquals(Object o)
	{
		Connection c = cast(Connection)o;
		return ((this.sourceNeuronIdx == c.sourceNeuronIdx)&&
				(this.targetNeuronIdx == c.targetNeuronIdx)&&
			 	(this.weight == c.weight)&&
				(this.signal == c.signal)&&
				(this.hebb_rate == c.hebb_rate)&&
				(this.hebb_pre_rate == c.hebb_pre_rate)
			);
	}
}
