<template>
    <div class="plotgrid-container">
        <div class="panel panel-default">
            <div class="panel-heading">
                <div class="row">
                    <div class="col-md-10 col-sm-9 col-xs-8">
                        <ul class="nav nav-pills">
                            <li>
                                <a href="#" class="msg-danger">错误
                                    <span class="badge badge-danger">42</span>
                                </a>
                            </li>
                            <li>
                                <a href="#" class="msg-warning">警告
                                    <span class="badge badge-warning">42</span>
                                </a>
                            </li>
                            <li>
                                <a href="#" class="msg-info">信息
                                    <span class="badge badge-info">3</span>
                                </a>
                            </li>
                        </ul>
                    </div>
                    <div class="col-md-2 col-sm-3 col-xs-4">
                        <div class="dropdown">
                            <button type="button" class="btn dropdown-toggle" id="dropdownMenu1" data-toggle="dropdown">选择时间
                                <span class="caret"></span>
                            </button>
                            <ul class="dropdown-menu" role="menu" aria-labelledby="dropdownMenu1">
                                <li role="presentation">
                                    <a role="menuitem" tabindex="-1" href="#">当前</a>
                                </li>
                                <li role="presentation">
                                    <a role="menuitem" tabindex="-1" href="#">一周</a>
                                </li>
                                <li role="presentation">
                                    <a role="menuitem" tabindex="-1" href="#">一月</a>
                                </li>
                                <li role="presentation">
                                    <a role="menuitem" tabindex="-1" href="#">三月</a>
                                </li>
                            </ul>
                        </div>
                    </div>
                </div>
            </div>
            <div class="panel-body">
                <div class="row">
                    <div class="col-md-8">
                        <h4>CPU监控</h4>
                        <div class="btn-group-vertical">
                            <button type='button' class='close-btn' data-dismiss='alert' title="移除">
                                <i class="fa fa-times-circle" aria-hidden="true"></i>
                            </button>
                            <button type="button" class="details-btn" title="查看详情">
                                <i class="fa fa-eye" aria-hidden="true"></i>
                            </button>
                        </div>
                        <img src="http://constey.de/wp-content/uploads/2013/06/esxi-nagios-cpu.png" class="img-thumbnail row-4">
                    </div>
                    <div class="col-md-4">
                        <h4>内存监控</h4>
                        <div class="btn-group-vertical">
                            <button type='button' class='close-btn' data-dismiss='alert'>
                                <i class="fa fa-times-circle" aria-hidden="true"></i>
                            </button>
                            <button type="button" class="details-btn">
                                <i class="fa fa-eye" aria-hidden="true"></i>
                            </button>
                        </div>
                        <img src="http://tse2.mm.bing.net/th?id=OIP.gywdsSJeYUfKh_WyxFC8MwEsBR&pid=15.1" class="img-thumbnail row-2">
                    </div>
                    <div class="col-md-4">
                        <h4>存储监控</h4>
                        <div class="btn-group-vertical">
                            <button type='button' class='close-btn' data-dismiss='alert'>
                                <i class="fa fa-times-circle" aria-hidden="true"></i>
                            </button>
                            <button type="button" class="details-btn">
                                <i class="fa fa-eye" aria-hidden="true"></i>
                            </button>
                        </div>
                        <img src="http://www.mfservis.cz/images/Nagios_graf_CPU.jpg" class="img-thumbnail row-2">            
                    </div>
                <!-- </div> -->
                <!-- <div class="row"> -->
                    <div class="col-md-4">
                        <h4>IB网络监控</h4>
                        <div class="btn-group-vertical">
                            <button type='button' class='close-btn' data-dismiss='alert'>
                                <i class="fa fa-times-circle" aria-hidden="true"></i>
                            </button>
                            <button type="button" class="details-btn">
                                <i class="fa fa-eye" aria-hidden="true"></i>
                            </button>
                        </div>
                        <img src="https://www.ibm.com/developerworks/cn/linux/l-ganglia-nagios-1/figure3.jpg" class="img-thumbnail row-4">
                    </div>
                    <div class="col-md-4">
                        <h4>网络监控</h4>
                        <div class="btn-group-vertical">
                            <button type='button' class='close-btn' data-dismiss='alert'>
                                <i class="fa fa-times-circle" aria-hidden="true"></i>
                            </button>
                            <button type="button" class="details-btn">
                                <i class="fa fa-eye" aria-hidden="true"></i>
                            </button>
                        </div>
                        <img src="http://filesimg.111cn.net/2014/04/17/20140417124113307.png" class="img-thumbnail row-2">
                    </div>
                    <div class="col-md-4">
                        <h4>电源监控</h4>
                        <div class="btn-group-vertical">
                            <button type='button' class='close-btn' data-dismiss='alert'>
                                <i class="fa fa-times-circle" aria-hidden="true"></i>
                            </button>
                            <button type="button" class="details-btn">
                                <i class="fa fa-eye" aria-hidden="true"></i>
                            </button>
                        </div>
                        <img src="http://tse2.mm.bing.net/th?id=OIP.HWWB8zLkO9PFPK-kL5Wi3gEsB_&pid=15.1" class="img-thumbnail row-2">            
                    </div>
                <!-- </div> -->
                <!-- <div class="row"> -->
                    <div class="col-md-4">
                        <div class="placeholder"></div>
                        <button type='button' onclick="showForm()" class='plus row-2' data-dismiss='alert'>+</button>
                    </div>
                </div>
            </div>
        </div>
    </div>
</template>

<script>
    // CSS Style
    import 'bootstrap/dist/css/bootstrap.css'
    import 'bootstrap/dist/css/bootstrap-theme.css'
    import 'font-awesome-webpack'
    export default {
        name: "plot-grid",
        data() {
            return {
                apiUrl: ""
            }
        }
    }
</script>

<style scoped>
    .msg-danger {
        background-color: #d9534f!important;
        color: #fff;
    }

    .msg-warning {
        background-color: #f0ad4e!important;
        color: #fff;
    }

    .msg-info {
        background-color: #337ab7!important;
        color: #fff;
    }
    .badge-danger {
        background-color: #fff!important;
        color: #000!important;           
    }
    .badge-warning {
        background-color: #fff!important;
        color: #000!important; 
    }
    .badge-info {
        background-color: #fff!important;
        color: #000!important;
    }
    .plotgrid-container {
        margin-top: 10px;
        width: 100%;
    }
    .panel {
        border: none;
        background-color: inherit;
    }

    .panel-heading {
        border-radius: 5px;
        background-color: #ddd!important;
    }
    .panel-body {
        margin-top: 10px;
        border: solid 1px #ddd;
        border-radius: 5px;
        background-color: #ddd;
    }
    .panel-body > .row {
        width: 99.99%;
        margin-left: -10px;
    }
    .btn {
        background-color: #fff;
    }
    .dropdown-toggle {
        padding: 9px 12px!important;
    }
    .img-thumbnail {
        width: 100%;
        height: 100%;
        margin-bottom: 10px;
    }
    .btn-group-vertical {
        position: absolute;
        right: -10px;
        opacity: 0.8;
        width: 26px;
        line-height: initial;
    }
    .close-btn {
        position: relative;
        padding: 2px 2px;
        float: left;
        width: 26px;
        display: block;
    }
    .details-btn {
        position: relative;
        padding: 2px 2px;
        float: left;
        width: 26px;
        display: block;
    }
    .plus {
        font-size: 60px;
        width: 100%;
        border-radius: 5px;
        border: solid 1px #ddd;
    }
    .row-2 {
        height: 200px;
    }
    .row-4 {
        height: 439px;
    }
    .row-6 {
        height: 639px;
    }
    h4 {
        text-align: center;
    }
    .placeholder {
        margin: 0px 10px 10px 10px;
        height: 19px;
    }
</style>