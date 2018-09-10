<template>
    <div class="container">
        <button type='button' class='close' data-dismiss='alert'>
            <i class="fa fa-times-circle" aria-hidden="true"></i>
        </button>
        <form :action="apiUrl" method="POST" role="form">
            <div class="form-group">
                <label for="monitor-index">选择指标</label>
                <select class="form-control" id="monitor-index">
                    <option>1</option>
                    <option>2</option>
                    <option>3</option>
                    <option>4</option>
                    <option>5</option>
                </select>
                <label for="visual-type">可视化类型</label>
                <select class="form-control" id="visual-type">
                    <option>1</option>
                    <option>2</option>
                    <option>3</option>
                    <option>4</option>
                    <option>5</option>
                </select>
            </div>
            <div class="form-group">
                <input type="text" class="form-control" placeholder="监控名称">
            </div>
            <div class="checkbox">
                <label><input type="checkbox" name="perm-storage">是否保存为配置?</label>
            </div>
            <button type="submit" class="btn btn-default">提交</button>
        </form>
    </div>
</template>

<script>
    // CSS Style
    import 'bootstrap/dist/css/bootstrap.css'
    import 'bootstrap/dist/css/bootstrap-theme.css'
    export default {
        name: "list-form",
        data() {
            return {
                apiUrl: ""
            }
        }
    }
</script>

<style scoped>
    .container {
        position: absolute;
        background-color: #ddd;
        top: 30%;
        left: 20%;
        width: 60%;
        padding: 10px;
        margin: auto auto;
        border-radius: 5px;
        border: solid 1px #fff;
    }    
</style>