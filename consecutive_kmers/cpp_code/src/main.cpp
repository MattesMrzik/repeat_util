#include "ConsecutiveKmers.h"
#include "argparser.h"

int main(int argc, char *argv[])
{
    Args args = parseArgs(argc, argv);
    ConsecutiveKmers ck(args);

    // #pragma omp parallel for
    for (auto &file_name : args.files)
    {
        ck.scan_file(file_name);
    }
    return 0;
}