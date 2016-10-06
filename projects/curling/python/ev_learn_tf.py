# -*- coding: utf-8 -*-
# ev_learn.py
# Katsuki Ohto
# http://kivantium.hateblo.jp/entry/2015/11/18/233834

import sys
import numpy as np

import tensorflow as tf

import dc
import image_log

#import model.linear as mdl
import model.cnn as mdl

in_size = dc.IMAGE_SIZE
out_size = dc.SCORE_LENGTH

N_PHASES = 16

if __name__ == "__main__":

    args = sys.argv
 
    # データのロードと設定
    image_logs = image_log.load_image_log(args[1])
    mepoch = int(args[2])
    batch_size = int(args[3])

    loaded_data_num = len(image_logs)

    print("number of loaded image-logs = %d" % loaded_data_num)
    
    # 場合分けのそれぞれのパターン数を求める
    # training data, test dataを作る
    phase_index_vector = [[] for _ in range(N_PHASES)]
    for i, il in enumerate(image_logs):
        #phase_index_vector[dc.to_turn_color(il[1])].append(i)
        phase_index_vector[il[1]].append(i)
    
    phase_data_num = [len(v) for v in phase_index_vector]
    data_num = np.sum(phase_data_num)

    print(phase_data_num)

    for ph in range(N_PHASES):
        np.random.shuffle(phase_index_vector[ph])
    
    image_vector = [np.empty((phase_data_num[ph], dc.IMAGE_WIDTH, dc.IMAGE_LENGTH, dc.IMAGE_PLAINS)
                             , dtype = np.float32) for ph in range(N_PHASES)]
    score_vector = [np.zeros((phase_data_num[ph], out_size), dtype = np.float32)
                    for ph in range(N_PHASES)]
    
    score_vector_sum = np.zeros(out_size, dtype = np.float32)

    for ph in range(N_PHASES):
        for i in range(phase_data_num[ph]):
            il = image_logs[phase_index_vector[ph][i]]
        
            end = il[0]
            turn = il[1]
        
            image_vector[ph][i] = il[4]
        
            sc_index = dc.StoIDX(il[3])
            score_vector[ph][i][sc_index] = 1
            score_vector_sum[sc_index] += 1 # 事前解析

    print("number of used image-logs = %d" % data_num)

    test_rate = 0.2
    test_num = []
    train_num = []
    
    for ph in range(N_PHASES):
        test_num.append(int(phase_data_num[ph] * test_rate))
        train_num.append(phase_data_num[ph] - test_num[ph])

    print("score distribution = ")
    print((score_vector_sum / data_num * 1000).astype(int))
        
    with tf.Graph().as_default():
        # 画像を入れる仮のTensor
        x_placeholder = tf.placeholder("float",
                                       shape = (None, dc.IMAGE_WIDTH, dc.IMAGE_LENGTH, dc.IMAGE_PLAINS))
        # ラベルを入れる仮のTensor
        labels_placeholder = tf.placeholder("float", shape = (None, dc.SCORE_LENGTH))
        # dropout率を入れる仮のTensor
        keep_prob = tf.placeholder("float")
        
        # モデル作成
        y_soft = [mdl.make(x_placeholder, keep_prob) for ph in range(N_PHASES)]
        # loss計算
        loss_value = [mdl.calc_loss(y_soft[ph], labels_placeholder) for ph in range(N_PHASES)]
        # training()を呼び出して訓練
        train_op = [mdl.train(loss_value[ph], 0.0005) for ph in range(N_PHASES)]
        # 精度の計算
        acc_op = [mdl.calc_accuracy(y_soft[ph], labels_placeholder) for ph in range(N_PHASES)]
        
        # 保存の準備
        saver = tf.train.Saver()
        # Sessionの作成
        sess = tf.Session()
        # 変数の初期化
        sess.run(tf.initialize_all_variables())
        # TensorBoardで表示する値の設定
        #summary_op = tf.merge_all_summaries()
        #summary_writer = tf.train.SummaryWriter(FLAGS.train_dir, sess.graph_def)
        
        for e in range(mepoch):
            
            # training phase
            for ph in range(N_PHASES):
                for bi in range((train_num[ph] - 1) // batch_size + 1):
                
                    # make batch
                    tmp_batch_size = min(batch_size, train_num[ph] - bi * batch_size)
                
                    istart = bi * batch_size
                    iend = min((bi + 1) * batch_size, train_num[ph])

                    # トレーニング
                    # feed_dictでplaceholderに入れるデータを指定する
                    sess.run(train_op[ph], feed_dict={
                         x_placeholder : image_vector[ph][istart : iend],
                         labels_placeholder : score_vector[ph][istart : iend],
                         keep_prob : 0.5})

                # test phase
                # 1 step終わるたびに精度を計算する
            for ph in range(N_PHASES):
                acc_sum = 0
                for i in range(train_num[ph], phase_data_num[ph], 1):
                    acc_sum += sess.run(acc_op[ph],
                                     feed_dict = {
                                      x_placeholder : image_vector[ph][i : (i + 1)],
                                      labels_placeholder : score_vector[ph][i : (i + 1)],
                                      keep_prob : 1.0})
                print("[epoch %d] phase %d, test accuracy %g" % (e, ph, acc_sum / test_num[ph]))
            
            # 1 step終わるたびにTensorBoardに表示する値を追加する
            """summary_str = sess.run(summary_op, feed_dict={
                                   x_placeholder : image_vector_all[],
                                   labels_placeholder: score_vector_all,
                                   keep_prob: 1.0})
                                   summary_writer.add_summary(summary_str, step)"""
            
            save_path = saver.save(sess, "ev_model_phase.ckpt")

    np.set_printoptions(threshold = 100000)
    #print(sess.run(W))
    #print(sess.run(b))


