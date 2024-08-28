⣿⣿⠛⠛⠛⠛⠿⠿⢿⠷⠸⠿⠿⠿⠿⢿⡿⠿⠛⠛⠛⠛⢻⣿⣿⣿⣿
⣿⣿⡷ ⣸⣟⠿⢿⢶⣶⣾⣿⣿⣿⣿⣷⣖⠻⢛⣿⢀⠰⣿⣿⣿⣿⣿
⣿⣿⣧⡐⢿⠟⠂⢈⣼⣿⣿⣿⣿⣿⣿⣿⣿⣖⠚⣍⠌⢠⣿⣿⣿⣿⣿
⣿⣿⣿⣿⣆⢰⣾⠉⢩⣽⠿⢿⣿⣿⣿⣿⠛⣤ ⡈⢰⣿⣿⣿⣿⣿⣿
⣿⣿⣿⣿⣿⣦⢻⣇⠘⠋ ⣘⣿⣿⣿⡇⠈⠉⢀⠇⣸⣿⣿⣿⣿⣿⣿
⣿⣿⣿⠟⠛⢃⣼⣿⣿⣷⣿⣿⣿⣿⣿⣿⣿⣿⠇⣠⣉⠛⠿⢿⣿⣿⣿
⣿⡿⠋⢠⡆⢸⣿⣿⣿⣿⣿⣿⡟⢛⠛⢻⣿⠋⢠⡟⢹⡶⢤⣌⣻⣿⡿
⠋⣤⣿⡿⢁⣾⣿⣿⣿⣿⣿⣿⣿⣟⣳⠞⣽⣄⢻⣽⠛⢠⣿⣿⣿⣿⣿
⣿⣿⡿⠁⣾⣿⣿⣿⣿⣿⣿⣿⣭⣯⣵⣿⣿⣿⡨⣿⣿⣿⣿⣿⣿⣿⣷
⣿⣿⡇⠰⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣧⢿⣿⣿⣿⣿⣿⣿⣿
⣿⣿⣷ ⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣼⣿⣿⣿⣿⣿⣿⣿
⣿⣿⣿⡃⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⢿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿

Simple Bose–Chaudhuri–Hocquenghem (BCH) Encoder / Decoder
So far supports creating a code that fixes up to 7 bit errors.
To fix more, lists have to be replaced with numpy arrays.
Can create galois fields. Tested up to GF(2^8). Fields of higher 
power might not work due to small ints in python. Numpy rewrite might be required.
