#!/bin/bash
# ./smiles2svg.py --smiles "c1ccccc1" --depiction 0 --hs > ./out/benzene_filled.svg
# ./smiles2svg.py --smiles "c1ccccc1" --depiction 1 --hs > ./out/benzene_ball_and_stick.svg
# ./smiles2svg.py --smiles "c1ccccc1" --depiction 2 --hs > ./out/benzene_wire.svg
# ./smiles2svg.py --smiles "c1ccccc1" --depiction 0 > ./out/benzene_filled_no_hs.svg
# ./smiles2svg.py --smiles "c1ccccc1" --depiction 1 > ./out/benzene_ball_and_stick_no_hs.svg
# ./smiles2svg.py --smiles "c1ccccc1" --depiction 2 > ./out/benzene_wire_no_hs.svg

# ./sdf2svg.py --sdf ./data/penicillin_G.sdf --depiction 0 --hs > ./out/penicillin_filled.svg
# ./sdf2svg.py --sdf ./data/penicillin_G.sdf --depiction 1 --hs > ./out/penicillin_ball_and_stick.svg
# ./sdf2svg.py --sdf ./data/penicillin_G.sdf --depiction 2 --hs > ./out/penicillin_wire.svg
# ./sdf2svg.py --sdf ./data/penicillin_G.sdf --depiction 0 > ./out/penicillin_filled_no_hs.svg
# ./sdf2svg.py --sdf ./data/penicillin_G.sdf --depiction 1 > ./out/penicillin_ball_and_stick_no_hs.svg
# ./sdf2svg.py --sdf ./data/penicillin_G.sdf --depiction 2 > ./out/penicillin_wire_no_hs.svg

./sdf2svg.py --sdf ./data/penicillin_G.sdf --depiction 0 --hs --yrot 0 > ./out/penicillin_filled_yrot0.svg
./sdf2svg.py --sdf ./data/penicillin_G.sdf --depiction 0 --hs --yrot 45 > ./out/penicillin_filled_yrot45.svg
./sdf2svg.py --sdf ./data/penicillin_G.sdf --depiction 0 --hs --yrot 90 > ./out/penicillin_filled_yrot90.svg
./sdf2svg.py --sdf ./data/penicillin_G.sdf --depiction 0 --hs --yrot 135 > ./out/penicillin_filled_yrot135.svg
./sdf2svg.py --sdf ./data/penicillin_G.sdf --depiction 0 --hs --yrot 180 > ./out/penicillin_filled_yrot180.svg
./sdf2svg.py --sdf ./data/penicillin_G.sdf --depiction 0 --hs --yrot 225 > ./out/penicillin_filled_yrot225.svg
./sdf2svg.py --sdf ./data/penicillin_G.sdf --depiction 0 --hs --yrot 270 > ./out/penicillin_filled_yrot270.svg
./sdf2svg.py --sdf ./data/penicillin_G.sdf --depiction 0 --hs --yrot 315 > ./out/penicillin_filled_yrot315.svg

./sdf2svg.py --sdf ./data/penicillin_G.sdf --depiction 1 --hs --yrot 0 > ./out/penicillin_ball_and_stick_yrot0.svg
./sdf2svg.py --sdf ./data/penicillin_G.sdf --depiction 1 --hs --yrot 45 > ./out/penicillin_ball_and_stick_yrot45.svg
./sdf2svg.py --sdf ./data/penicillin_G.sdf --depiction 1 --hs --yrot 90 > ./out/penicillin_ball_and_stick_yrot90.svg
./sdf2svg.py --sdf ./data/penicillin_G.sdf --depiction 1 --hs --yrot 135 > ./out/penicillin_ball_and_stick_yrot135.svg
./sdf2svg.py --sdf ./data/penicillin_G.sdf --depiction 1 --hs --yrot 180 > ./out/penicillin_ball_and_stick_yrot180.svg
./sdf2svg.py --sdf ./data/penicillin_G.sdf --depiction 1 --hs --yrot 225 > ./out/penicillin_ball_and_stick_yrot225.svg
./sdf2svg.py --sdf ./data/penicillin_G.sdf --depiction 1 --hs --yrot 270 > ./out/penicillin_ball_and_stick_yrot270.svg
./sdf2svg.py --sdf ./data/penicillin_G.sdf --depiction 1 --hs --yrot 315 > ./out/penicillin_ball_and_stick_yrot315.svg

./sdf2svg.py --sdf ./data/penicillin_G.sdf --depiction 2 --hs --yrot 0 > ./out/penicillin_wire_yrot0.svg
./sdf2svg.py --sdf ./data/penicillin_G.sdf --depiction 2 --hs --yrot 45 > ./out/penicillin_wire_yrot45.svg
./sdf2svg.py --sdf ./data/penicillin_G.sdf --depiction 2 --hs --yrot 90 > ./out/penicillin_wire_yrot90.svg
./sdf2svg.py --sdf ./data/penicillin_G.sdf --depiction 2 --hs --yrot 135 > ./out/penicillin_wire_yrot135.svg
./sdf2svg.py --sdf ./data/penicillin_G.sdf --depiction 2 --hs --yrot 180 > ./out/penicillin_wire_yrot180.svg
./sdf2svg.py --sdf ./data/penicillin_G.sdf --depiction 2 --hs --yrot 225 > ./out/penicillin_wire_yrot225.svg
./sdf2svg.py --sdf ./data/penicillin_G.sdf --depiction 2 --hs --yrot 270 > ./out/penicillin_wire_yrot270.svg
./sdf2svg.py --sdf ./data/penicillin_G.sdf --depiction 2 --hs --yrot 315 > ./out/penicillin_wire_yrot315.svg

./sdf2svg.py --sdf ./data/penicillin_G.sdf --depiction 0 > ./out/penicillin_filled_no_hs.svg
./sdf2svg.py --sdf ./data/penicillin_G.sdf --depiction 1 > ./out/penicillin_ball_and_stick_no_hs.svg
./sdf2svg.py --sdf ./data/penicillin_G.sdf --depiction 2 > ./out/penicillin_wire_no_hs.svg

./color_substructure.py --sdf ./data/penicillin_G.sdf --depiction 0 --hs > ./out/penicillin_G_filled_lactam.svg
./color_substructure.py --sdf ./data/penicillin_G.sdf --depiction 1 --hs > ./out/penicillin_G_ball_and_stick_lactam.svg
./color_substructure.py --sdf ./data/penicillin_G.sdf --depiction 2 --hs > ./out/penicillin_G_wire_lactam.svg

./color_substructure.py --sdf ./data/penicillin_G.sdf --depiction 0 > ./out/penicillin_G_filled_crippen_logp.svg
./color_substructure.py --sdf ./data/penicillin_G.sdf --depiction 1 > ./out/penicillin_G_ball_and_stick_crippen_logp.svg
./color_substructure.py --sdf ./data/penicillin_G.sdf --depiction 2 > ./out/penicillin_G_wire_crippen_logp.svg