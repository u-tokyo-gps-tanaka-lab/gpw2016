mkdir -p ./record/mf_p/
../../projects/curling/out/release/server_gat -g 800 -e 2 -t 80000000 -p 9876 -l ./record/mf_p/ -s 1 2>./record/mf_p/output_server &
../../projects/curling/out/release/client_mcts_full -npo 200000 -p 9876 -id mcts_full -s 1 2>./record/mf_p/output_mcts_full &
../../projects/curling/out/release/client_pure -npo 200000 -p 9876 -id pure -s 1 2>./record/mf_p/output_pure &
mkdir -p ./record/mpf_pp/
../../projects/curling/out/release/server_gat -g 800 -e 2 -t 80000000 -p 9877 -l ./record/mpf_pp/ -s 2 2>./record/mpf_pp/output_server &
../../projects/curling/out/release/client_mcts_policy_full -npo 200000 -p 9877 -id mcts_policy_full -s 2 2>./record/mpf_pp/output_mcts_policy_full &
../../projects/curling/out/release/client_pure_policy -npo 200000 -p 9877 -id pure_policy -s 2 2>./record/mpf_pp/output_pure_policy &
mkdir -p ./record/mpf_mf/
../../projects/curling/out/release/server_gat -g 800 -e 2 -t 80000000 -p 9888 -l ./record/mpf_mf/ -s 3 2>./record/mpf_mf/output_server &
../../projects/curling/out/release/client_mcts_policy_full -npo 200000 -p 9888 -id mcts_policy_full -s 3 2>./record/mpf_mf/output_mcts_policy_full &
../../projects/curling/out/release/client_mcts_full -npo 200000 -p 9888 -id mcts_full -s 3 2>./record/mpf_mf/output_mcts_full &
mkdir -p ./record/pp_p/
../../projects/curling/out/release/server_gat -g 800 -e 2 -t 80000000 -p 9879 -l ./record/pp_p/ -s 4 2>./record/pp_p/output_server &
../../projects/curling/out/release/client_pure_policy -npo 200000 -p 9879 -id pure_policy -s 4 2>./record/pp_p/output_pure_policy &
../../projects/curling/out/release/client_pure -npo 200000 -p 9879 -id pure -s 4 2>./record/pp_p/output_pure &