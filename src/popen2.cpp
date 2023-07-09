#include "popen2.h"
#include <unistd.h>
#include <iostream>

pid_t
popen2(const char *command, const char* pars, int *infp, int *outfp)
{
    int p_stdin[2], p_stdout[2];
    pid_t pid;

    if (pipe(p_stdin) != 0 || pipe(p_stdout) != 0)
        return -1;

    pid = fork();

    if (pid < 0)
        return pid;
    else if (pid == 0)
    {
        close(p_stdin[WRITE]);
        dup2(p_stdin[READ], READ);
        close(p_stdout[READ]);
        dup2(p_stdout[WRITE], WRITE);

        //execl("/bin/sh", "sh" "-c", command, NULL);
        //execl("/bin/sh", "sh" "-c", //"../htd --child-limit 1", NULL);
	//const char* addr[] = {command, "--child-limit", "1", "--opt", "width", "--iterations", "50"};
        execl(
		command, command, /*const_cast<char* const *>(addr)); //command,*/ "--child-limit", "2", 
		"--opt", "width", "--iterations", "50000", "--patience", "1500", "--strategy", "min-fill",  NULL);
	//std::cerr << command << std::endl;
        perror("execl");
        exit(1);
    }

    if (infp == NULL)
        close(p_stdin[WRITE]);
    else
        *infp = p_stdin[WRITE];

    if (outfp == NULL)
        close(p_stdout[READ]);
    else
        *outfp = p_stdout[READ];

    return pid;
}
