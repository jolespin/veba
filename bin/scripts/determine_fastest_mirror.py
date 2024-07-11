#!/usr/bin/env python
import sys, os, argparse, requests, time
# from tqdm import tqdm

__program__ = os.path.split(sys.argv[0])[-1]
__version__ = "2024.6.5"

# # URLs to test
# url_mirror_1="https://data.gtdb.ecogenomic.org/releases/release214/214.1/auxillary_files/gtdbtk_r214_data.tar.gz"
# url_mirror_2="https://data.ace.uq.edu.au/public/gtdb/data/releases/release214/214.1/auxillary_files/gtdbtk_r214_data.tar.gz"

def download_and_measure_speed(url, test_duration):
    start_time = time.time()
    end_time = start_time + test_duration
    total_bytes = 0

    with requests.get(url, stream=True) as r:
        r.raise_for_status()
        for chunk in r.iter_content(chunk_size=8192):
            if time.time() >= end_time:
                break
            total_bytes += len(chunk)
    
    elapsed_time = time.time() - start_time
    speed_bps = total_bytes * 8 / elapsed_time  # Convert bytes to bits
    speed_kbps = speed_bps / 1024  # Convert bits per second to kilobits per second
    return speed_kbps
    
def main(args=None):
    # Path info
    script_directory  =  os.path.dirname(os.path.abspath( __file__ ))
    script_filename = __program__

    # Path info
    description = """
    Running: {} v{} via Python v{} | {}""".format(__program__, __version__, sys.version.split(" ")[0], sys.executable)
    usage = "{} -t 10 url_1 url_2".format(__program__)
    epilog = "Copyright 2021 Josh L. Espinoza (jespinoz@jcvi.org)"

    # Parser
    parser = argparse.ArgumentParser(description=description, usage=usage, epilog=epilog, formatter_class=argparse.RawTextHelpFormatter)

    # Pipeline
    parser.add_argument("urls", nargs="+", help = "URLs to check speed. Must include at least 1.")
    parser.add_argument("-t","--timeout", type=int, default=10, help = "Number of seconds to use for each speed test [Default: 10]")
    parser.add_argument("-m","--mode", type=str, default="fastest", choices={"table", "fastest"}, help = "Output format [Default: fastest]")


    # Options
    opts = parser.parse_args()
    opts.script_directory  = script_directory
    opts.script_filename = script_filename
    
    # Test all URLs
    url_to_speed = dict()
    for url in opts.urls:
        url_to_speed[url] = download_and_measure_speed(url, opts.timeout)
        
    # Write URLs to stderr
    if opts.mode == "table":
        f_table = sys.stdout
    if opts.mode == "fastest":
        f_table = sys.stderr
    print("URL", "Speed [kb/s]", sep="\t", file=f_table)
    for url in opts.urls:
        print(url, url_to_speed[url], sep="\t", file=f_table)
        
    # Print fastest
    if opts.mode == "fastest":
        print(sorted(url_to_speed.items(), key=lambda x: x[1], reverse=False)[0][0], file=sys.stdout)
        


if __name__ == "__main__":
    main()
    
                

