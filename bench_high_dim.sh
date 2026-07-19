#! /bin/bash

metrics=(0 1 2)
ns=(12 16 18)
dims=(16 32 64)
deltas=(16 64 256 1024)

printf "[Prot] [Metric] [Dim] [Delta] [Size] [Com.(MB)] [Time(s)]
"

for metric in "${metrics[@]}"; do
  for n in "${ns[@]}"; do
    for dim in "${dims[@]}"; do
      for delta in "${deltas[@]}"; do
        if [[ "$metric" == "2" && "$n" == "18" && "$dim" == "64" ]]; then
          continue
        fi
        ./build/fpsi -high -d $dim -delta $delta -nn $n -p $metric -try 1
      done
    done
  done
done
