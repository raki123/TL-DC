// TLDC "Too long; Didn't Count" A length limited path counter.
// Copyright (C) 2023 Rafael Kiesel, Markus Hecher

// This program is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.

// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.

// You should have received a copy of the GNU General Public License
// along with this program.  If not, see <https://www.gnu.org/licenses/>.

#include "popen2.h"
#include <iostream>
#include <unistd.h>

pid_t popen2(const char *command, const char *pars, int *infp, int *outfp) {
  int p_stdin[2], p_stdout[2];
  pid_t pid;

  if (pipe(p_stdin) != 0 || pipe(p_stdout) != 0)
    return -1;

  pid = fork();

  if (pid < 0)
    return pid;
  else if (pid == 0) {
    close(p_stdin[WRITE]);
    dup2(p_stdin[READ], READ);
    close(p_stdout[READ]);
    dup2(p_stdout[WRITE], WRITE);

    // execl("/bin/sh", "sh" "-c", command, NULL);
    // execl("/bin/sh", "sh" "-c", //"../htd --child-limit 1", NULL);
    // const char* addr[] = {command, "--child-limit", "1", "--opt", "width",
    // "--iterations", "50"};
    execl(command, command,
          /*const_cast<char* const *>(addr)); //command,*/ "--child-limit", "2",
          "--opt", "width", "--iterations", "50000", "--patience", "1500",
          "--strategy", "min-fill", NULL);
    // std::cerr << command << std::endl;
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
