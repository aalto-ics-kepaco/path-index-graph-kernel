package mapper;

import mapper.AstarMapping;
import java.util.*;

class FeatDiffComparator implements Comparator<AstarMapping>
{
	public int compare(AstarMapping a, AstarMapping b)
	{
		if (a.getFeatDiff() < b.getFeatDiff())
			return -1;
		else if (a.getFeatDiff() > b.getFeatDiff())
			return 1;
		return 0;
	}
}

