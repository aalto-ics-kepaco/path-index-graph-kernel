package mechanism;

public class Timer
{
	private long startTime = 0;
	private long stopTime = 0;

	public Timer()
	{
	}

	public Timer(boolean init)
	{
		if (init)
			this.startTime = System.currentTimeMillis();
	}

	public void start()
	{
		this.startTime = System.currentTimeMillis();
	}

	public long getTime()
	{
		if (this.stopTime != 0)
			return this.stopTime - this.startTime;
		if (this.startTime != 0)
			return System.currentTimeMillis() - this.startTime;
		else
			return -1;
	}

	public void stop()
	{
		this.stopTime = System.currentTimeMillis();
	}

	public void reset()
	{
		this.startTime = 0;
		this.stopTime = 0;
	}
}

