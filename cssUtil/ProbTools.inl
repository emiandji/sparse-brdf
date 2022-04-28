



template <typename TypeSig, typename TypeProb>
void HistCount(const std::vector<TypeSig>& vSig, const std::vector<TypeProb>& vEdges, std::vector<TypeProb>& vProb)
{
	if (vProb.empty()) //add exception here
		return;
	std::vector<size_t> vn(vProb.size(), 0);

	size_t sum = 0;
	for (size_t i = 0; i < vSig.size(); i++)
	{
		for (size_t j = 0; j < vEdges.size() - 1; j++)
		{
			if ((TypeProb)vSig[i] > vEdges[j] && (TypeProb)vSig[i] < vEdges[j + 1])
			{
				vn[j]++;
				sum++;
			}
		}
	}

	for (size_t i = 0; i < vProb.size(); i++)
		vProb[i] = (TypeProb)vn[i] / (TypeProb)sum;
}