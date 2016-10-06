# -*- coding: utf-8 -*-
# cnn.py
# Katsuki Ohto
# https://www.tensorflow.org/versions/0.6.0/tutorials/mnist/pros/index.html
# http://kivantium.hateblo.jp/entry/2015/11/18/233834

import numpy as np
import tensorflow as tf

import dc

N_FILTERS = 64

# convolutional neural network
def make(x_placeholder, keep_prob):
    # モデルを作成する関数

    # 初期化関数
    def weight_variable(shape):
        initial = tf.truncated_normal(shape, stddev = 0.01)
        return tf.Variable(initial)

    def bias_variable(shape):
        initial = tf.constant(0.0, shape = shape)
        return tf.Variable(initial)

    # 畳み込み層
    def conv2d(x, W):
        return tf.nn.conv2d(x, W, strides = [1, 1, 1, 1], padding = 'SAME')

    # プーリング層
    def max_pool_2x2(x):
        return tf.nn.max_pool(x, ksize = [1, 2, 2, 1],
                          strides = [1, 2, 2, 1], padding = 'SAME')
    # inputの変形
    # x_image = tf.reshape(x_placeholder, [-1, 28, 28, 1])

    # 畳み込み層1
    with tf.name_scope('conv1') as scope:
        W_conv1 = weight_variable([5, 5, 1, N_FILTERS])
        b_conv1 = bias_variable([N_FILTERS])
        print(tf.shape(W_conv1))
        print(tf.shape(x_placeholder))
        print(tf.shape(b_conv1))
        h_conv1 = tf.nn.relu(conv2d(x_placeholder, W_conv1) + b_conv1)

    # 畳み込み層2
    with tf.name_scope('conv2') as scope:
        W_conv2 = weight_variable([3, 3, N_FILTERS, N_FILTERS])
        b_conv2 = bias_variable([N_FILTERS])
        h_conv2 = tf.nn.relu(conv2d(h_conv1, W_conv2) + b_conv2)

    # 畳み込み層3
    with tf.name_scope('conv3') as scope:
        W_conv3 = weight_variable([3, 3, N_FILTERS, N_FILTERS])
        b_conv3 = bias_variable([N_FILTERS])
        h_conv3 = tf.nn.relu(conv2d(h_conv2, W_conv3) + b_conv3)

    # 畳み込み層4
    with tf.name_scope('conv4') as scope:
        W_conv4 = weight_variable([3, 3, N_FILTERS, N_FILTERS])
        b_conv4 = bias_variable([N_FILTERS])
        h_conv4 = tf.nn.relu(conv2d(h_conv3, W_conv4) + b_conv4)

    # 全結合層1
    with tf.name_scope('fc1') as scope:
        W_fc1 = weight_variable([28 * 28 * N_FILTERS, 256])
        b_fc1 = bias_variable([256])
        h_pool3_flat = tf.reshape(h_conv4, [-1, 28 * 28 * N_FILTERS])
        h_fc1 = tf.nn.relu(tf.matmul(h_pool3_flat, W_fc1) + b_fc1)
        # dropout
        h_fc1_drop = tf.nn.dropout(h_fc1, keep_prob)

    # 全結合層2
    with tf.name_scope('fc2') as scope:
        W_fc2 = weight_variable([256, dc.SCORE_LENGTH])
        b_fc2 = bias_variable([dc.SCORE_LENGTH])

    # softmax
    with tf.name_scope('softmax') as scope:
        y = tf.nn.softmax(tf.matmul(h_fc1_drop, W_fc2) + b_fc2)

    return y # 正規化されいていない確率配列を返す

def normalize(y_soft):
    return y_soft / tf.reduce_sum(y_soft)

def calc_loss(y_soft, labels):
    cross_entropy = -tf.reduce_sum(labels * tf.log(tf.clip_by_value(y_soft, 1e-10, 1.0)))
    # TensorBoardで表示するよう指定
    #tf.scalar_summary("cross_entropy", cross_entropy)
    return cross_entropy

def train(loss, learning_rate = 0.0001):
    train_step = tf.train.AdamOptimizer(learning_rate).minimize(loss)
    #train_step = tf.train.GradientDescentOptimizer(1.0).minimize(cross_entropy)
    return train_step

def calc_accuracy(y_soft, labels):
    correct_prediction = tf.equal(tf.argmax(y_soft, 1), tf.argmax(labels, 1))
    accuracy = tf.reduce_mean(tf.cast(correct_prediction, "float"))
    #tf.scalar_summary("accuracy", accuracy)
    return accuracy
