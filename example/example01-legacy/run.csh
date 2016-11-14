${SHELL_DIR}/bin/h2cp-legacy /home/mcaprio/research/data/h2/run0164-ob-9/JISP16-ob-9-20.bin tbme-VNN.bin 2 2
${SHELL_DIR}/bin/h2cp-legacy /home/mcaprio/research/data/h2/run0164-ob-9/VC-ob-9-20.bin tbme-VC.bin 2 2

${SHELL_DIR}/bin/h2gen-legacy < h2gen-legacy.in
${SHELL_DIR}/bin/h2mixer-legacy < h2mixer-legacy.in
