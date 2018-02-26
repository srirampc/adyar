# adyar
Approximate Common Substring

To compile program:
`
g++ -std=c++11 -O3 -DNDEBUG -I ~/include -L ~/lib kkmacs_v2.cpp -o kkmacsv2 -lsdsl -ldivsufsort -ldivsufsort64
`

To run:
`./kkmacsv2 <text_file_name> <value_of_k>`

Example:
`./kkmacsv2 text.txt 1`

Output:
Runtime, value of symmetric acs distance
