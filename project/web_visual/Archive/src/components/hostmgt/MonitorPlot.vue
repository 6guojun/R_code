<template>
    <div class="container">
        <plot-grid></plot-grid>
        <list-form v-if="isActive"></list-form>
    </div>
</template>

<script>
    import PlotGrid from "../common/PlotGrid"
    import ListForm from "../common/ListForm"

    export default {
        name: 'monitor-plot',
        props: {
            apiUrl: "",
            isActive: false
        },
        components: {
            PlotGrid,
            ListForm
        }
    }
</script>

<style scoped>
    .container {
        width: 100%;
    }
</style>