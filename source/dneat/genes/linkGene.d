module dneat.genes.linkGene;

public class LinkGene
{
private:

	// The IDs of the neurons that this link connects
	uint _fromNeuronID, _toNeuronID; 
   	// The link's innovation ID
	uint _innovationID;
	// This variable is modified during evolution
	// The weight of the connection
	double _weight;

	// Is it recurrent?
	bool _isRecurrent;

public:

	// Weight of the connection
	double weight() const @property
	{
		return _weight;
	}
	void weight(double weight) @property
	{
		_weight = weight;
	}	

	this(uint inID, uint outID, uint innovID, double wgt, bool reccurent = false)
	{
		_fromNeuronID = inID;
		_toNeuronID = outID;
		_innovationID = innovID;
		_weight = wgt;
		_isRecurrent = reccurent;
	}
	this(){}
	
	uint fromNeuronID() const @property
	{
		return _fromNeuronID;
	}

	uint toNeuronID() const @property
	{
		return _toNeuronID;
	}

	uint innovationID() const @property
	{
		return _innovationID;
	}

	bool isRecurrent() const @property
	{
		return _isRecurrent;
	}

	bool isLoopedRecurrent() const @property
	{
		return _fromNeuronID == _toNeuronID;
	}

	// Overloading Operators using Innovation ID as Criteria

	override bool opEquals(Object g)
	{
		auto lg = cast(LinkGene)g;
		return _innovationID == lg._innovationID;
	}

	override int opCmp(Object g)
	{
		auto lg = cast(LinkGene)g;
		if(_innovationID < lg._innovationID)
		{
			return -1;
		}
		else if(_innovationID > lg._innovationID)
		{
			return 1;
		}
		else
		{
			return 0;
		}
	}
}
