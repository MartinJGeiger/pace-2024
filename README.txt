This is our contribution to the 'Heuristic'-track of the Parameterized Algorithms 
and Computational Experiments Challenge PACE 2024.

The program has been written in C# programming langue using  Microsoft Visual 
Studio 2019, complied with both VS and Mono 6.12.0.206, and tested under Windows
and Linux. 
At this point, we believe that the code works, at least for the test instances
provided by the challenge organizers.



compile it as such:
mcs -platform:x64 -optimize+ -out:pace-2024.exe pace-2024.cs



There is a single important parameter, the 'MaxRuntime', which is set to 300. 
As C# does not handle SIGTERM signals, the program keeps checking it's runtime, 
and terminates before (MaxRuntime - 1) seconds. 
This behaviour has been verified in more than 24,000 runs on an Intel(R) Xeon(R) 
Platinum 8268 CPU @ 2.90GHz (hence a CPU not too different to the Intel(R) 
Xeon(R) Gold 6342 CPU @ 2.80GHz). 
Note that on optil.io, in order to avoid a ""Time Limit Exceeded"-error (TLE), the
MaxRuntime must be set to around 294.


Martin J. Geiger, m.j.geiger@hsu-hh.de, JUNE 2024