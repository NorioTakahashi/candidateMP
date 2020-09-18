## NT4 の変更の覚書

### Sept 18, 2020 by norio



* NT1のctrlファイルの読み込みは

!!ad_comm::change_datafile_name("NT1_ver0_ctrl.dat"); // set parameter values for MP

でやっているので（注意!!!!：このファイル名はmycontrol*.datのMP実行を指定する行の
ファイル名と一致していないといけない）

NT2でctrlを**"NT2_ctrl.dat"**として読み込むためには

  !!ad_comm::change_datafile_name(**"NT2_ctrl.dat"**); // set parameter values for MP

とソースコードを変更しておかないといけない。

* ↑であったのでNT4でもそうなっていたが、この入力ファイル読み込みの行をコメントアウトすることでmycontrol*.datのMP実行を指定する行で、ファイル名も任意に指定できるようにすることができることが分かったので変更した。

  **//**  !!ad_comm::change_datafile_name("NT4_ctrl.dat"); // set parameter values for MP

* 変更は**NT4.tpl**と**NT4_2020.tpl**（2020年のCTPによるプロジェクション結果との比較のために使ったNT4）に加えた。